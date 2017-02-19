import os
import numpy as np
import nipy as ni
import nipype as npe
import nibabel as nib
from glob import glob

from nipype.interfaces import fsl
from nipype.interfaces import spm

def normalize(brain_data):
    """
  	Uses SPM to normalize the volume to 2mm cubic voxel.
    """
	pass

def coregister(ref_volume,target_volume):
    """
  	Uses SPM to coregister the two inputed volumes.
    """	
    pass

def segment_anat(anat_volume):
    """
  	Uses SPM to segment GM, WM, and CSF tissue.
    """	
    pass

def register_to_standard(brain_data, standard_brain):
    """
  	Affine registration to a template brain.
    """
	pass

def mask_funct(brain_volume,mask_volume):
    """
  	Uses fslmaths to remove all non-GM voxels.
    """
    #fslmaths wrcbf.nii -mas c1wspgr.nii wrcbf_gm.nii
	pass

