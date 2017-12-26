
# coding: utf-8

# In[ ]:

from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.io import FreeSurferSource, SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node

from nipype.interfaces.freesurfer import Binarize, MRIConvert, FSCommand
from nipype.interfaces.fsl import ApplyMask, Reorient2Std
from nipype.interfaces.fsl.preprocess import MCFLIRT, SliceTimer, FLIRT, FAST
from nipype.interfaces.spm.preprocess import Smooth
from nipype.algorithms.rapidart import ArtifactDetect
from nipype.interfaces.fsl.model import GLM
from nipype.algorithms.confounds import CompCor
from nipype.interfaces.nipy.preprocess import Trim

# MATLAB setup - Specify path to current SPM and the MATLAB's default mode
from nipype.interfaces.matlab import MatlabCommand
MatlabCommand.set_default_paths('~/spm12/toolbox')
MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")

#other study-specific variables
project_home = '/Users/catcamacho/Dropbox/Projects/TH_NAR_ASL/proc'
raw_dir = project_home + '/raw'
subjects_list = open(project_home + '/misc/subjects.txt').read().splitlines()
output_dir = project_home + '/proc/preprocessing'
wkflow_dir = project_home + '/workflows'
template = project_home + '/template/MNI152_T1_2mm_brain.nii'

#freesurfer setup
subjects_dir = project_home + '/freesurfer'
FSCommand.set_default_subjects_dir(subjects_dir)

#Population specific variables for ASL
nex_asl = 3 #number of excitations from the 3D ASL scan parameters
inversion_efficiency = 0.8 #from GE
background_supp_eff = 0.75 #from GE
efficiency = inversion_efficiency * background_supp_eff 
T1_blood = 1.6 #T1 of blood in seconds(1.6s at 3T and 1.4s at 1.5T)
sat_time = 2 #in seconds, from GE
partition_coeff = 0.9 #whole brain average in ml/g
scaling_factor = 32 #scaling factor, can be taken from PW dicom header at position 0043,107f (corresponds to #coils?)
postlabel_delay = 1.525 #post label delay in seconds
labeling_time = 1.450 #labeling time in seconds
T1_tissue = 1.2 #estimated T1 of grey matter in seconds
TR = 4.844 #repetition time

# Study variables for resting state processing
rest_TR = 2 #in seconds
num_slices = 29
slice_direction = 3 #3 = z direction
interleaved = True
#all rates are in Hz (1/TR or samples/second)
highpass_freq = 0.008 #in Hz
lowpass_freq = 0.09 #in Hz
vols_to_trim = 6

smoothing_kernel = 6 #in mm


# In[ ]:

## File handling nodes

# Select subjects
infosource = Node(IdentityInterface(fields=['subjid']),
                  name='infosource')
infosource.iterables = [('subjid', subjects_list)]

# SelectFiles
templates = {'pw_volume': raw_dir + '/{subjid}/pw.nii',
             'rest': raw_dir + '/{subjid}/rest_raw.nii',
             'pd_volume': raw_dir + '/{subjid}/pd.nii'}
selectfiles = Node(SelectFiles(templates), name='selectfiles')

# FreeSurferSource - Data grabber specific for FreeSurfer data
fssource = Node(FreeSurferSource(subjects_dir=subjects_dir),
                run_without_submitting=True,
                name='fssource')
# Datasink
datasink = Node(DataSink(base_directory = output_dir, 
                         container = output_dir), 
                name='datasink')

# DataSink output substitutions (for ease of folder naming)
substitutions = [('_subjid_', '')]
datasink.inputs.substitutions = substitutions

# Grab processed anat
struct_template = {'proc_anat':output_dir + '/preproc_anat/{subjid}/brainmask_out_reoriented.nii'}
structgrabber = Node(SelectFiles(struct_template),name='structgrabber')


# In[ ]:

## File Processing nodes

# convert files to nifti
mri_convert = Node(MRIConvert(out_type='nii',
                              conform_size=2,
                              crop_size= (128, 128, 128)), 
                   name='mri_convert')

# reorient data for consistency
reorient_anat = Node(Reorient2Std(output_type='NIFTI'),
                     name='reorient_anat')

# Binarize -  binarize and dilate image to create a brainmask
binarize = Node(Binarize(min=0.5,
                         max=300,
                         dilate=3,
                         erode=2,
                         out_type='nii'),
                name='binarize')

# reorient PD/PW data for consistency
reorient_pd = Node(Reorient2Std(output_type='NIFTI'),
                   name='reorient_pd')

reorient_pw = Node(Reorient2Std(output_type='NIFTI'),
                   name='reorient_pw')

# Register subject's PD/PW to subject anatomy
reg_pw2anat = Node(FLIRT(output_type='NIFTI',
                         out_file='warped_pw.nii'),
                   name='reg_pw2anat')

reg_pd2anat = Node(FLIRT(output_type='NIFTI',
                         apply_xfm = True, 
                         out_file='warped_pd.nii'),
                   name='reg_pd2anat')

# Mask brain in pw volume
applyMask_pw = Node(ApplyMask(output_type='NIFTI'),
                    name='applyMask_pw')

# Mask brain pd volumes
applyMask_pd = Node(ApplyMask(output_type='NIFTI'), 
                    name='applyMask_pd')


# In[ ]:

## Custom functions
#quantify CBF from PW volume (Alsop MRIM 2015 method)
def quantify_cbf_alsop(pw_volume,pd_volume,sat_time,postlabel_delay,T1_blood,labeling_time,efficiency,partition_coeff,TR,T1_tissue,scaling_factor,nex_asl):
    import os
    import nibabel as nib
    from numpy import exp
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    
    # set variables
    pw_nifti1 = nib.load(pw_volume)
    pw_data = pw_nifti1.get_data()
    pw_data = pw_data.astype(float)
    pd_nifti1 = nib.load(pd_volume)
    pd_data = pd_nifti1.get_data()
    pd_data = pd_data.astype(float)
    conversion = 6000 #to convert values from mL/g/s to mL/100g/min
    pd_factor = 1/(1-exp((-1*TR)/T1_tissue))
    
    cbf_numerator = conversion*partition_coeff*pw_data*exp(postlabel_delay/T1_blood)
    cbf_denominator = sat_time*efficiency*T1_blood*scaling_factor*nex_asl*pd_data*pd_factor*(1-exp((-1*labeling_time)/T1_blood))
    cbf_data = cbf_numerator/cbf_denominator
    
    cbf_volume = nib.Nifti1Image(cbf_data, pw_nifti1.affine)
    nib.save(cbf_volume, 'alsop_cbf.nii')
    cbf_path = os.path.abspath('alsop_cbf.nii')
    return cbf_path
    
def make_histogram(cbf_vol):
    from os.path import abspath
    from nibabel import load
    from numpy import isnan
    from matplotlib import pyplot as plt
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    
    cbf_image = load(cbf_vol)
    cbf_data = cbf_image.get_data()

    plt.Figure()
    data = cbf_data.flatten()
    plt.hist(data[~isnan(data)], bins = range(0,200,10))
    plt.title("Distribution of CBF values in mL/100g/min")
    plt.savefig('cbf_histogram.png')
    cbf_hist = abspath('cbf_histogram.png')
    return(cbf_hist)

def create_cbfcoreg_plot(ref,anat):
    from os.path import abspath
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    from nilearn import plotting
    
    coreg_filename='coregistration.png'
    display = plotting.plot_anat(ref, display_mode='ortho',
                                 draw_cross=False,
                                 title = 'coregistration to anatomy')
    display.add_edges(anat)
    display.savefig(coreg_filename) 
    display.close()
    coreg_file = abspath(coreg_filename)
    
    return(coreg_file)


# In[ ]:

## Normalizing data for first and second level analysis

quant_cbf_alsop = Node(Function(input_names=['pw_volume','pd_volume',
                                             'sat_time','postlabel_delay',
                                             'T1_blood','labeling_time',
                                             'efficiency','partition_coeff',
                                             'TR','T1_tissue','scaling_factor',
                                             'nex_asl'],
                                output_names=['cbf_volume'],
                                function=quantify_cbf_alsop),
                       name='quant_cbf_alsop')
quant_cbf_alsop.inputs.sat_time=sat_time
quant_cbf_alsop.inputs.postlabel_delay=postlabel_delay
quant_cbf_alsop.inputs.T1_blood=T1_blood
quant_cbf_alsop.inputs.labeling_time=labeling_time
quant_cbf_alsop.inputs.efficiency=efficiency
quant_cbf_alsop.inputs.partition_coeff=partition_coeff
quant_cbf_alsop.inputs.TR=TR
quant_cbf_alsop.inputs.T1_tissue=T1_tissue
quant_cbf_alsop.inputs.scaling_factor=scaling_factor
quant_cbf_alsop.inputs.nex_asl=nex_asl

# make histgram of CBF values for each subject
makecbf_hist = Node(Function(input_names=['cbf_vol'], 
                             output_names=['cbf_hist'], 
                             function=make_histogram), 
                    name = 'makecbf_hist')

checkcoreg_cbf = Node(Function(input_names=['ref','anat'],
                               output_names=['coreg_file'],
                               function=create_cbfcoreg_plot), 
                      name='checkcoreg_cbf')

checkcoreg_mni = Node(Function(input_names=['ref','anat'],
                               output_names=['coreg_file'],
                               function=create_cbfcoreg_plot), 
                      name='checkcoreg_mni')
checkcoreg_mni.inputs.anat=template

# Register subject's anatomy/CBF to the template
reg_anat2mni = Node(FLIRT(output_type='NIFTI',
                          reference=template),
                    name='reg_anat2mni')

reg_cbf2mni = Node(FLIRT(output_type='NIFTI',
                         apply_xfm = True,
                         reference=template,
                         out_file='warped_cbf.nii'),
                   name='reg_cbf2mni')

# smooth cbf data
cbf_smooth = Node(Smooth(fwhm=[6,6,6], 
                     implicit_masking=True), 
              name='smooth')


# In[ ]:

# Create a flow for preprocessing anat + asl volumes 
cbfpreproc = Workflow(name='cbfpreproc')

# Connect all components of the preprocessing workflow
cbfpreproc.connect([(infosource,selectfiles, [('subjid', 'subjid')]),
                    (infosource,fssource, [('subjid','subject_id')]),
                    (fssource,mri_convert, [('brainmask', 'in_file')]),
                    (mri_convert,reorient_anat, [('out_file','in_file')]),
                    (reorient_anat,binarize, [('out_file','in_file')]),
                    (reorient_anat,reg_pw2anat, [('out_file','reference')]),
                    (reorient_anat,reg_pd2anat, [('out_file','reference')]),
                    (reorient_anat,checkcoreg_cbf, [('out_file','anat')]),
                    (reorient_anat,reg_anat2mni, [('out_file','in_file')]),
                    (reg_anat2mni,reg_cbf2mni, [('out_matrix_file','in_matrix_file')]),
                    (reg_anat2mni,checkcoreg_mni, [('out_file','ref')]),
                    (binarize,applyMask_pw, [('binary_file','mask_file')]),
                    (binarize,applyMask_pd, [('binary_file','mask_file')]),
                    
                    (selectfiles,reorient_pw, [('pw_volume','in_file')]),
                    (reorient_pw,reg_pw2anat, [('out_file','in_file')]),
                    (reg_pw2anat,reg_pd2anat, [('out_matrix_file','in_matrix_file')]),
                    (selectfiles,reorient_pd, [('pd_volume','in_file')]),
                    (reorient_pd,reg_pd2anat, [('out_file','in_file')]),
                    (reg_pw2anat,applyMask_pw, [('out_file','in_file')]),
                    (applyMask_pw, quant_cbf_alsop, [('out_file','pw_volume')]),
                    (reg_pd2anat,applyMask_pd, [('out_file','in_file')]),
                    (applyMask_pd,quant_cbf_alsop, [('out_file','pd_volume')]),
                    
                    (quant_cbf_alsop,checkcoreg_cbf, [('cbf_volume','ref')]),
                    (quant_cbf_alsop,reg_cbf2mni,[('cbf_volume','in_file')]),
                    (reg_cbf2mni, cbf_smooth, [('out_file','in_files')]),
                    
                    (checkcoreg_cbf,datasink, [('coreg_file','CBF_anat_coreg')]),
                    (checkcoreg_mni,datasink, [('coreg_file','MNI_anat_coreg')]),
                    (cbf_smooth,datasink, [('smoothed_files','std_cbf')]),
                    (reorient_anat,datasink, [('out_file','preproc_anat')]),
                    (reg_anat2mni,datasink, [('out_file','std_anat')]),
                    (reg_anat2mni,datasink, [('out_matrix_file','reg_mni_xform')]),
                    (binarize,datasink, [('binary_file','brainmask')])                    
                   ])
cbfpreproc.base_dir = wkflow_dir
cbfpreproc.write_graph(graph2use='flat')
cbfpreproc.run('MultiProc', plugin_args={'n_procs': 2})
