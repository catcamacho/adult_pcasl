
# coding: utf-8

# In[2]:

import os
from glob import glob
import nibabel as nib
from numpy import exp
from os.path import join as opj

from nipype.interfaces import freesurfer as fs
from nipype.interfaces.ants import Registration, ApplyTransforms
from nipype.interfaces.c3 import C3dAffineTool
from nipype.interfaces import fsl
from nipype.interfaces.utility import Function, IdentityInterface, Merge
from nipype.interfaces.io import FreeSurferSource, SelectFiles, DataSink
from nipype.interfaces.slicer.registration import brainsresample
from nipype.algorithms.misc import Gunzip
from nipype.pipeline.engine import Workflow, Node, MapNode

# MATLAB setup - Specify path to current SPM and the MATLAB's default mode
from nipype.interfaces.matlab import MatlabCommand
MatlabCommand.set_default_paths('~/spm12/toolbox')
MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")

#other study-specific variables
project_home = '/Users/catcamacho/Dropbox/Projects/TH_NAR_ASL/proc/dev'
subjects_dir = project_home
raw_dir = project_home + '/raw'
subjects_list = os.listdir(raw_dir)
output_dir = project_home + '/proc'
template = project_home + '/template/MNI152_T1_2mm_brain.nii'
gray_matter_mask = project_home + '/template/MNI152_T1_2mm_graymatter.nii'

#freesurfer setup
fs_dir = subjects_dir + '/freesurfer'
fs.FSCommand.set_default_subjects_dir(fs_dir)

#Population specific variables for ASL
nex_asl = 3 #number of excitations from the 3D ASL scan parameters
inversion_efficiency = 0.8 #from GE
background_supp_eff = 0.75 #from GE
efficiency = inversion_efficiency * background_supp_eff 
T1_blood = 1.6 #T1 of blood in seconds(1.6s at 3T and 1.4s at 1.5T)
sat_time = 2 #in seconds, from GE
partition_coeff = 0.9 #whole brain average in ml/g
scaling_factor = 32 #scaling factor, can be taken from PW dicom header at position 0043,107f
postlabel_delay = 2.025 #1.525 #post label delay in seconds
labeling_time = 1.450 #labeling time in seconds
T1_tissue = 1.2 #estimated T1 of grey matter in seconds


# In[3]:

## File handling nodes

# Select subjects
infosource = Node(IdentityInterface(fields=['subjid', 'volume']),
                  name='infosource')
infosource.iterables = [('subjid', subjects_list),('volume',['pw','pd'])]


# SelectFiles
templates = {'asl_volume': raw_dir + '/{subjid}/{volume}.nii.gz'}
selectfiles = Node(SelectFiles(templates), name='selectfiles')

# FreeSurferSource - Data grabber specific for FreeSurfer data
fssource = Node(FreeSurferSource(subjects_dir=fs_dir),
                run_without_submitting=True,
                name='fssource')
# Datasink
datasink = Node(DataSink(), name='datasink')
datasink.inputs.base_directory = output_dir
datasink.inputs.container = output_dir
# DataSink output substitutions (for ease of folder naming)
substitutions = [('_subjid_', ''),
                ('volume_',''),
                ('_reoriented',''),
                ('_warped','')]
datasink.inputs.substitutions = substitutions


# In[4]:

## File Processing nodes

# convert files to nifti
mri_convert = Node(fs.MRIConvert(out_type='nii',
                                out_orientation='RAS',
                                conform_size=2,
                                crop_size= (128, 128, 128)), 
                   name='mri_convert')

# reorient data for consistency
reorient = Node(fsl.utils.Reorient2Std(output_type='NIFTI'),
                name='reorient')

# BBRegister - coregister a volume to the Freesurfer anatomical
bbregister = Node(fs.BBRegister(init='header',
                                contrast_type='t2',
                                out_fsl_file=True),
                  name='bbregister')

# Volume Transformation - transform the brainmask into functional space
applyVolTrans = Node(fs.ApplyVolTransform(inverse=True),
                     name='applyVolTrans')

# Binarize -  binarize and dilate image to create a brainmask
binarize = Node(fs.Binarize(min=0.5,
                         dilate=2,
                         out_type='nii'),
                name='binarize')

# Mask brain in pw and pd volumes
applyMask = Node(fsl.maths.ApplyMask(output_type='NIFTI'), 
                        name='applyMask')

# N3 bias correction using MINC tools, will be necessary for babies
nu_correct = Node(fs.MNIBiasCorrection(), name='nu_correct')


# In[11]:

# Create a flow for preprocessing anat + asl volumes 
preprocflow = Workflow(name='preprocflow')

# Connect all components of the preprocessing workflow
preprocflow.connect([(infosource, selectfiles, [('subjid', 'subjid')]),
                     (infosource, selectfiles, [('volume', 'volume')]),
                     (infosource, fssource, [('subjid','subject_id')]),
                     (infosource, bbregister, [('subjid','subject_id')]),
                     (fssource, mri_convert, [('brainmask', 'in_file')]),
                     (mri_convert, datasink, [('out_file', 'brainmask_nifti')]),
                     (mri_convert, binarize, [('out_file', 'in_file')]),
                     (binarize, datasink, [('binary_file', 'binary_mask')]),
                     (binarize, applyMask, [('binary_file','mask_file')]),
                     (selectfiles, reorient, [('asl_volume', 'in_file')]),
                     (reorient, bbregister, [('out_file', 'source_file')]),
                     (infosource, applyVolTrans, [('subjid','subject')]),
                     (bbregister, applyVolTrans, [('out_reg_file', 'reg_file')]),
                     (reorient, applyVolTrans, [('out_file', 'target_file')]),
                     (mri_convert, applyVolTrans, [('out_file', 'source_file')]),
                     (applyVolTrans, datasink, [('transformed_file', 'registered_volumes')]),
                     (applyVolTrans, applyMask, [('transformed_file','in_file')]),
                     (applyMask, datasink, [('out_file','masked_volumes')])
                    ])
preprocflow.base_dir = opj(output_dir)
preprocflow.write_graph(graph2use='flat')
preprocflow.run('MultiProc', plugin_args={'n_procs': 2})


# In[6]:

## File handling nodes for CBF proc

# Select subjects
cbfinfosource = Node(IdentityInterface(fields=['subjid']),
                  name='cbfinfosource')
cbfinfosource.iterables = [('subjid', subjects_list)]


# SelectFiles
templates = {'pw_volume': output_dir + '/masked_volumes/{subjid}_pw/pw_masked.nii',
            'pd_volume': output_dir + '/masked_volumes/{subjid}_pd/pd_masked.nii',
            'anat_volume': output_dir + '/brainmask_nifti/{subjid}_pw/brainmask_out.nii'}
cbfselectfiles = Node(SelectFiles(templates), name='cbfselectfiles')

# Datasink
cbfdatasink = Node(DataSink(), name='cbfdatasink')
cbfdatasink.inputs.base_directory = output_dir
cbfdatasink.inputs.container = output_dir
# DataSink output substitutions (for ease of folder naming)
substitutions = [('_subjid_', '')]
cbfdatasink.inputs.substitutions = substitutions


# In[7]:

## Custom functions

#quantify CBF from PW volume
def quantify_cbf(pw_volume,pd_volume,sat_time,T1_tissue,postlabel_delay,T1_blood,labeling_time,efficiency,nex_asl,scaling_factor,partition_coeff):
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
    
    cbf_numerator = (1-exp(-1*sat_time/T1_tissue))*exp(postlabel_delay/T1_blood)
    cbf_denominator = 2*T1_blood*(1-exp(-1*labeling_time/T1_blood))*efficiency*nex_asl
    cbf_data = conversion*partition_coeff*(cbf_numerator/cbf_denominator)*(pw_data/(scaling_factor*pd_data))
    
    cbf_volume = nib.Nifti1Image(cbf_data, pw_nifti1.affine)
    nib.save(cbf_volume, 'cbf.nii')
    cbf_path = os.path.abspath('cbf.nii')
    return cbf_path
   
    

quant_cbf = Node(name='quant_cbf',
                interface=Function(input_names=['pw_volume','pd_volume',
                                                'sat_time','T1_tissue',
                                                'postlabel_delay','T1_blood',
                                                'labeling_time','efficiency',
                                                'nex_asl','scaling_factor','partition_coeff'],
                                  output_names=['cbf_volume'],
                                  function=quantify_cbf))
quant_cbf.inputs.sat_time=sat_time
quant_cbf.inputs.T1_tissue=T1_tissue
quant_cbf.inputs.postlabel_delay=postlabel_delay
quant_cbf.inputs.T1_blood=T1_blood
quant_cbf.inputs.labeling_time=labeling_time
quant_cbf.inputs.efficiency=efficiency
quant_cbf.inputs.nex_asl=nex_asl
quant_cbf.inputs.scaling_factor=scaling_factor
quant_cbf.inputs.partition_coeff=partition_coeff


# In[8]:

## Normalizing data for first and second level analysis

# Register subject's anatomy to the template
linearReg = Node(fsl.FLIRT(output_type='NIFTI',
                          reference=template),
                         name='linearReg')

# Robust registration- did worse than FNIRT.
robustReg = Node(fs.RobustRegister(target_file=template,
                                  terminal_output='file',
                                  init_orient=True,
                                  auto_sens=True),
                name='robustReg')

applyxform = Node(fs.ApplyVolTransform(inverse=False,
                                      target_file=template),
                     name='applyxform')

## Register CBF vol to MNI space
# Volume Transformation - transform the cbf volume into MNI space
warpCBF = Node(fs.ApplyVolTransform(inverse=False,
                                   target_file=template),
                     name='warpCBF')

# Mask Gray Matter
maskGrayMatter = Node(fsl.maths.ApplyMask(mask_file=gray_matter_mask,
                                         output_type='NIFTI'),
                     name='maskGrayMatter')


# In[9]:

# create a flow for quantifying CBF and warping to MNI space.
cbfprocflow = Workflow(name='cbfprocflow')

# connect the nodes
cbfprocflow.connect([(cbfinfosource, cbfselectfiles, [('subjid', 'subjid')]),
                     (cbfselectfiles, quant_cbf, [('pw_volume', 'pw_volume')]),
                     (cbfselectfiles, quant_cbf, [('pd_volume', 'pd_volume')]),
                     (quant_cbf, cbfdatasink, [('cbf_volume', 'quantified_cbf')]),
                     (cbfselectfiles, linearReg, [('anat_volume', 'in_file')]),
                     (linearReg, cbfdatasink, [('out_file','linwarped_anat')]),
                     (linearReg, warpCBF, [('out_matrix_file', 'fsl_reg_file')]),
                     (quant_cbf, warpCBF, [('cbf_volume', 'source_file')]),
                     (warpCBF, cbfdatasink, [('transformed_file', 'warped_cbf_vol')]),
                     (warpCBF, maskGrayMatter, [('transformed_file', 'in_file')]),
                     (maskGrayMatter, cbfdatasink, [('out_file', 'masked_cbf')])
                    ]),
cbfprocflow.base_dir = opj(output_dir)
cbfprocflow.write_graph(graph2use='flat')
cbfprocflow.run('MultiProc', plugin_args={'n_procs': 2})


# In[10]:

#### Not being used till ITK bug worked out! ####
# Register subject's anatomy to the template (more accurate version, to get transform for cbf)
antsreg = Node(Registration(args='--float',
                            collapse_output_transforms=True,
                            fixed_image=template,
                            initial_moving_transform_com=True,
                            num_threads=1,
                            output_inverse_warped_image=True,
                            output_warped_image=True,
                            sigma_units=['vox']*3,
                            transforms=['Rigid', 'Affine', 'SyN'],
                            terminal_output='file',
                            winsorize_lower_quantile=0.005,
                            winsorize_upper_quantile=0.995,
                            convergence_threshold=[1e-06],
                            convergence_window_size=[10],
                            metric=['MI', 'MI', 'CC'],
                            metric_weight=[1.0]*3,
                            number_of_iterations=[[1000, 500, 250, 100],
                                                  [1000, 500, 250, 100],
                                                  [100, 70, 50, 20]],
                            radius_or_number_of_bins=[32, 32, 4],
                            sampling_percentage=[0.25, 0.25, 1],
                            sampling_strategy=['Regular',
                                               'Regular',
                                               'None'],
                            shrink_factors=[[8, 4, 2, 1]]*3,
                            smoothing_sigmas=[[3, 2, 1, 0]]*3,
                            transform_parameters=[(0.1,),
                                                  (0.1,),
                                                  (0.1, 3.0, 0.0)],
                            use_histogram_matching=True,
                            write_composite_transform=True),
               name='antsreg')

## Register CBF vol to MNI space

# First get transform matrix:
# Coregister the median to the surface
bbreg = Node(fs.BBRegister(init='fsl',
                             contrast_type='t2',
                             out_fsl_file=True),
                  name='bbreg')

# Convert the BBRegister transformation to ANTS ITK format
convert2itk = Node(C3dAffineTool(fsl2ras=True,
                                 itk_transform=True),
                   name='convert2itk')


# Concatenate BBRegister's and ANTS' transforms into a list
merge = Node(Merge(2), iterfield=['in2'], name='mergexfm')

# Then do the warp
warpCBF = Node(ApplyTransforms(args='--float',
                                input_image_type=3,
                                interpolation='Linear',
                                invert_transform_flags=[False, False],
                                num_threads=1,
                                reference_image=template,
                                terminal_output='file'),
                name='warpCBF')

