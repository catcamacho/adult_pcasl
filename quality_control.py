import os
import numpy as np
import nipy as ni
import nipype as npe
import nibabel as nib
from glob import glob

def plot_hist(brain_volume):
    """
  	Plots a histogram of the inputed volume's voxel values.
    """
    pass
	
def save_slices_overlay(anat_volume,overlay_volume):
    """
  	Saves a slice by slice view of the overlay volume overlaid on the anat brain.
  	This is helpful for checking subject-specific masking, group-wise T-maps, etc.
    """
	pass
	
def plot_groups(group_list,volume_file):
    """
  	Creates a scatterplot of all mean whole-brain CBFs to visually inspect outliers.
    """
	pass