
# coding: utf-8

# In[ ]:


from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.io import FreeSurferSource, SelectFiles, DataSink
from nipype.pipeline.engine import Workflow, Node

from nipype.interfaces.freesurfer import Binarize, MRIConvert, FSCommand
from nipype.interfaces.fsl import ApplyMask, Reorient2Std, MotionOutliers
from nipype.interfaces.fsl.preprocess import MCFLIRT, SliceTimer, FLIRT, FAST, SUSAN
from nipype.algorithms.rapidart import ArtifactDetect
from nipype.interfaces.fsl.model import GLM
from nipype.algorithms.confounds import CompCor
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.nipy.preprocess import Trim


# MATLAB setup - Specify path to current SPM and the MATLAB's default mode
from nipype.interfaces.matlab import MatlabCommand
MatlabCommand.set_default_paths('Users/catcamacho/spm12/toolbox')
MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")

#other study-specific variables
project_home = '/Users/catcamacho/Dropbox/projects/th_nar_asl/proc'
raw_dir = project_home + '/raw'
output_dir = project_home + '/proc/rest_preproc'
workflow_dir = project_home + '/workflows'
asl_dir = project_home + '/proc/asl_preproc'
template = project_home + '/template/MNI152_T1_2mm_brain.nii'

subjects_list = open(project_home + '/misc/subjects.txt').read().splitlines()
#subjects_list = ['003-DT2']

#freesurfer setup
subjects_dir = project_home + '/freesurfer'
FSCommand.set_default_subjects_dir(subjects_dir)

# FSL set up- change default file output type
from nipype.interfaces.fsl import FSLCommand
FSLCommand.set_default_output_type('NIFTI_GZ')

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


# Select subjects
infosource = Node(IdentityInterface(fields=['subjid']),
                  name='infosource')
infosource.iterables = [('subjid', subjects_list)]

# SelectFiles
templates = {'rest': raw_dir + '/{subjid}/rest_raw.nii.gz', 
             'proc_anat':asl_dir + '/preproc_anat/{subjid}/brainmask_out_reoriented.nii',
             'mask':asl_dir + '/brainmask/{subjid}/brainmask_out_reoriented_thresh.nii'}
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


# In[ ]:


def bandpass_filter(in_file, lowpass, highpass, TR):
    import numpy as np
    import nibabel as nb
    from os.path import abspath
    from nipype.interfaces.afni.preprocess import Bandpass
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    
    out_file = 'func_filtered'
    bp = Bandpass()
    bp.inputs.highpass = highpass
    bp.inputs.lowpass = lowpass
    bp.inputs.in_file = in_file
    bp.inputs.tr = TR
    bp.inputs.out_file = out_file
    bp.inputs.outputtype = 'NIFTI'
    bp.run()
    
    out_file = abspath(out_file + '.nii.gz')
    return(out_file)

def adjust_masks(masks):
    from os.path import abspath
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    
    from nipype.interfaces.freesurfer.model import Binarize
    #pve0 = csf, pve1 = gm, pve2 = wm
    
    origvols = sorted(masks)
    csf = origvols[0]
    wm = origvols[2]
    vols = []
    
    wm_file = 'WM_seg.nii.gz'
    binary = Binarize()
    binary.inputs.in_file = wm
    binary.inputs.min = 0.5
    binary.inputs.max = 2
    binary.inputs.binary_file = wm_file
    binary.run()
    wm_new = abspath(wm_file)
    vols.append(wm_new)
    
    csf_file = 'CSF_seg.nii.gz'
    binary2 = Binarize()
    binary2.inputs.in_file = csf
    binary2.erode = 1
    binary2.inputs.min = 0.5
    binary2.inputs.max = 2
    binary2.inputs.binary_file = csf_file
    binary2.run()
    csf_new = abspath(csf_file)
    vols.append(csf_new)
    
    return(vols)

## This reads in the output from fsl_motion_outliers as opposed to ART
def create_noise_matrix(vols_to_censor,motion_params,comp_noise):
    from numpy import genfromtxt, zeros, column_stack, savetxt
    from os import path
    
    motion = genfromtxt(motion_params, delimiter=None, dtype=None, skip_header=0)
    comp_noise = genfromtxt(comp_noise, delimiter=None, dtype=None, skip_header=1)
    censor_vol_list = genfromtxt(vols_to_censor, delimiter=None, dtype=None, skip_header=0)
    
    (a,b) = censor_vol_list.shape
    if b > 0:
        scrubbing = censor_vol_list
        noise_matrix = column_stack((motion,comp_noise,scrubbing))
    else:
        try:
            c = censor_vol_list.size
        except:
            c = 0
        d=len(comp_noise)

        if c > 1:
            scrubbing = zeros((d,c),dtype=int)
            for t in range(c):
                scrubbing[censor_vol_list[t],t] = 1
            noise_matrix = column_stack((motion,comp_noise,scrubbing))
        elif c == 1:
            scrubbing = zeros((d,c),dtype=int)
            scrubbing[censor_vol_list] = 1
            noise_matrix = column_stack((motion,comp_noise,scrubbing))
        else:
            noise_matrix = column_stack((motion,comp_noise))
    
    noise_file = 'noise_matrix.txt'
    savetxt(noise_file, noise_matrix)
    noise_filepath = path.abspath(noise_file)
    
    return(noise_filepath)


# In[ ]:


def create_coreg_plot(epi,anat):
    from os.path import abspath
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    from nilearn import plotting
    from nipype.interfaces.nipy.preprocess import Trim
    
    epiVol = 'firstVol.nii.gz'
    trim = Trim()
    trim.inputs.in_file = epi
    trim.inputs.out_file = epiVol
    trim.inputs.end_index = 1
    trim.inputs.begin_index = 0
    trim.run()
    
    coreg_filename='coregistration.png'
    display = plotting.plot_anat(epiVol, display_mode='ortho',
                                 draw_cross=False,
                                 title = 'coregistration to anatomy')
    display.add_edges(anat)
    display.savefig(coreg_filename) 
    display.close()
    coreg_file = abspath(coreg_filename)
    
    return(coreg_file)

def check_mask_coverage(epi,brainmask):
    from os.path import abspath
    from nipype import config, logging
    config.enable_debug_mode()
    logging.update_logging(config)
    from nilearn import plotting
    from nipype.interfaces.nipy.preprocess import Trim
    
    epiVol = 'firstVol.nii.gz'
    trim = Trim()
    trim.inputs.in_file = epi
    trim.inputs.out_file = epiVol
    trim.inputs.end_index = 1
    trim.inputs.begin_index = 0
    trim.run()
    
    maskcheck_filename='maskcheck.png'
    display = plotting.plot_anat(epiVol, display_mode='ortho',
                                 draw_cross=False,
                                 title = 'check brainmask coverage')
    display.add_contours(brainmask,levels=[.5], colors='r')
    display.savefig(maskcheck_filename) 
    display.close()
    checkmask_file = abspath(maskcheck_filename)
    
    return(checkmask_file)
    
make_coreg_img = Node(Function(input_names=['epi','anat'],
                                         output_names=['coreg_file'],
                                         function=create_coreg_plot),
                      name='make_coreg_img')

make_checkmask_img = Node(Function(input_names=['epi','brainmask'],
                                         output_names=['maskcheck_file'],
                                         function=check_mask_coverage),
                          name='make_checkmask_img')


# In[ ]:


# Resting state preprocessing nodes

# Trim the initial volumes
trim = Node(Trim(begin_index=vols_to_trim), name='trim')

# reorient the functional data
reorient_func = Node(Reorient2Std(), name='reorient_func')

# Realign each volume to first volume: in_file; out_file, par_file
realign = Node(MCFLIRT(out_file='realigned_func.nii.gz',
                       save_plots=True, 
                       mean_vol=True
                      ), 
               name='realign')

# Slice time correction: in_file, slice_time_corrected_file
slicetime = Node(SliceTimer(time_repetition=rest_TR, 
                            interleaved=interleaved,
                            slice_direction=slice_direction, 
                            out_file='stfunc.nii.gz'), 
                    name='slicetime')

# register the functional volumes to the subject space anat
# inputs: in_file, reference; out_file out_matrix_file
reg_func_to_anat = Node(FLIRT(out_matrix_file='xform.mat'),
                        name='reg_func_to_anat')

apply_reg_to_func = Node(FLIRT(apply_xfm=True, 
                               out_file='warped_func.nii.gz'),  
                         name='apply_reg_to_func')

# Apply binary mask to merged functional scan: in_file, mask_file; out_file
mask_func = Node(ApplyMask(out_file='masked_func.nii.gz'), 
                 name='mask_func')

# Bandpass Filtering all rates are in Hz (1/TR or samples/second)
bandpass = Node(name='bandpass', 
                interface=Function(input_names=['in_file','lowpass','highpass','TR'], 
                                   output_names=['out_file'],
                                   function=bandpass_filter))
bandpass.inputs.lowpass = lowpass_freq
bandpass.inputs.highpass = highpass_freq
bandpass.inputs.TR = rest_TR

# gunzip the fmri file before putiing it through ART
unzip = Node(Gunzip(), name='unzip')

# Artifact detection for scrubbing/motion assessment
#art = Node(ArtifactDetect(mask_type='file',
#                          parameter_source='FSL',
#                          norm_threshold=0.25, #mutually exclusive with rotation and translation thresh
#                          zintensity_threshold=3,
#                          use_differences=[True, False]),
#           name='art')

get_FD = Node(MotionOutliers(threshold=0.2, 
                             metric='fd', 
                             no_motion_correction=False, 
                             out_file='outliers.txt', 
                             out_metric_plot='fd.png', 
                             out_metric_values='fd.txt'), 
              name='get_FD')

# Segment structural scan
segment = Node(FAST(no_bias=True, 
                    segments=True, 
                    number_classes=3), 
               name='segment')

# Fix the segmentations
fix_confs = Node(name='fix_confs',
                 interface=Function(input_names=['masks'], 
                                    output_names=['vols'],
                                    function=adjust_masks))
# actually run compcor
compcor = Node(CompCor(merge_method='none', 
                       num_components=3), 
               name='compcor')

# Create a denoising mask with compcor + motion
noise_mat = Node(name='noise_mat', interface=Function(input_names=['vols_to_censor','motion_params','comp_noise'],
                                                      output_names=['noise_filepath'], 
                                                      function=create_noise_matrix))

# Denoise the data
denoise = Node(GLM(out_res_name='denoised_residuals.nii.gz', 
                   out_data_name='denoised_func.nii.gz'), 
               name='denoise')


# In[ ]:


# workflow
## need to add in FD rather than ART

restpreproc = Workflow(name='restpreproc')
restpreproc.connect([(infosource,selectfiles,[('subjid','subjid')]),
                     (selectfiles,reg_func_to_anat,[('proc_anat','reference')]),
                     (selectfiles,apply_reg_to_func,[('proc_anat','reference')]),
                     (selectfiles,mask_func,[('mask','mask_file')]),
                     (selectfiles,segment, [('proc_anat','in_files')]),
                     
                     (selectfiles,reorient_func,[('rest','in_file')]),
                     (reorient_func,trim,[('out_file','in_file')]),
                     (trim,realign,[('out_file','in_file')]),
                     (realign, slicetime,[('out_file','in_file')]),
                     (slicetime,reg_func_to_anat,[('slice_time_corrected_file','in_file')]),
                     (slicetime,apply_reg_to_func,[('slice_time_corrected_file','in_file')]),
                     (reg_func_to_anat,apply_reg_to_func,[('out_matrix_file','in_matrix_file')]),
                     (apply_reg_to_func,mask_func,[('out_file','in_file')]),
                     
                     (segment,fix_confs,[('tissue_class_files','masks')]),
                     (fix_confs,compcor,[('vols','mask_files')]),
                     (trim,get_FD, [('out_file','in_file')]),
                     (get_FD, noise_mat,[('out_file','vols_to_censor')]),
                     (get_FD, noise_mat,[('out_metric_values','motion_params')]),
                     (mask_func,compcor,[('out_file','realigned_file')]),
                     (compcor,noise_mat,[('components_file','comp_noise')]),
                     (noise_mat,denoise,[('noise_filepath','design')]),
                     (mask_func,denoise,[('out_file','in_file')]),
                     (denoise,bandpass,[('out_data','in_file')]),
                     
                     (mask_func,make_coreg_img,[('out_file','epi')]),
                     (selectfiles,make_coreg_img,[('proc_anat','anat')]),
                     
                     (make_coreg_img,datasink,[('coreg_file','coregcheck_image')]),
                     (get_FD, datasink, [('out_metric_plot','FD_plot')]),
                     (bandpass,datasink,[('out_file','preproc_func')])        
                    ])
restpreproc.base_dir = workflow_dir
restpreproc.write_graph(graph2use='flat')
restpreproc.run('MultiProc', plugin_args={'n_procs': 2})

