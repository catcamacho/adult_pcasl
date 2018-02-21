
# coding: utf-8

# In[ ]:


import os
from os.path import join as opj
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink, DataGrabber
from nipype.interfaces.fsl.model import Randomise, Cluster
from nipype.interfaces.fsl.utils import Merge, ImageMeants
from nipype.interfaces.fsl.maths import ApplyMask
from nipype.interfaces.freesurfer.model import Binarize

#other study-specific variables
project_home = '/Users/catcamacho/Dropbox/Projects/TH_NAR_ASL/proc'
preproc_dir = project_home + '/proc/preprocessing'
output_dir = project_home + '/proc/secondlevel'
wkflow_dir = project_home + '/workflows'
mask = project_home + '/template/MNI_2mm_GM_mask.nii'

#covariate_file = project_home + '/modelinfo/MCageCov.txt'
t_contrasts = project_home + '/misc/tcon.con'
group_mat = project_home + '/misc/design.mat'


# In[ ]:


# Data handling nodes
grabcbfdata = Node(DataGrabber(template=preproc_dir + '/std_cbf/*/swarped_cbf.nii', 
                               sort_filelist=True, 
                               outfields=['cbf_list']), 
                   name='grabcbf')

# Datasink
datasink = Node(DataSink(base_directory = output_dir, 
                         container = output_dir), 
                name='datasink')

# DataSink output substitutions (for ease of folder naming)
substitutions = [('_subjid_', '')]
datasink.inputs.substitutions = substitutions


# In[ ]:


merge_cbf = Node(Merge(dimension = 't'), name='merge_cbf')

# FSL randomise for higher level analysis
highermodel = Node(Randomise(tfce=False,
                             c_thresh=2,
                             tcon=t_contrasts,
                             raw_stats_imgs= True,
                             mask=mask,
                             num_perm= 5000,
                             design_mat=group_mat),
                   name = 'highermodel')

# Make mask of significant voxels
binarize = MapNode(Binarize(min=0.95), name='binarize', iterfield=['in_file'])

# Mask the tstat files
applyMask = MapNode(ApplyMask(), name='applyMask',iterfield=['mask_file','in_file'])

cluster = MapNode(Cluster(threshold=2, 
                          out_index_file=True, 
                          out_localmax_txt_file=True, 
                          minclustersize=True, 
                          peak_distance=6), 
                  name='cluster', 
                  iterfield=['in_file'])


# In[ ]:


cbf_groupflow = Workflow(name='cbf_groupflow')
cbf_groupflow.connect([(grabcbfdata,merge_cbf, [('cbf_list','in_files')]),
                       (merge_cbf,highermodel, [('merged_file','in_file')]),
                       (highermodel, binarize, [('t_corrected_p_files','in_file')]),
                       (binarize, applyMask, [('binary_file','mask_file')]),
                       (highermodel, applyMask, [('tstat_files','in_file')]),
                       (applyMask, cluster, [('out_file','in_file')]),
                       (cluster, datasink, [('index_file','clusters'), 
                                            ('localmax_txt_file','localmax_txt_file')]),
                       (highermodel,datasink, [('t_corrected_p_files','t_corrected_p_files')]),
                       (highermodel,datasink, [('tstat_files','tstat_files')]),
                       (applyMask, datasink, [('out_file','masked_tstats')])
                      ])
cbf_groupflow.base_dir = wkflow_dir
cbf_groupflow.write_graph(graph2use='flat')
cbf_groupflow.run('MultiProc', plugin_args={'n_procs': 2})

