import os
import numpy as np
import nibabel as nib
import nipy as nip
import nipype as npe
from nipype.interfaces.nipy.preprocess import Trim
trim = Trim()

path = '/Users/catcamacho/Box Sync/CBFBIRN_practice/210-NART2'
imgFile = path + '/raw/asl.nii.gz'
img = nib.load(imgFile)

# Determine file type and transform to NIFTI if need be

# Determine file properties
def getImgMeta(img):
	header = img.header
	dim = header.get_data_shape()
	matrix = dim[0:2]
	numSlices = dim[2]
	numVols = dim[3]
	
	return matrix
	return numSlices
	return numVols

# Separate the asl and pd volumes
def splitASLvols(imgFile, path):
	aslVolFile = path + '/aslVol.nii.gz'
	pdVolFile = path + '/pdVol.nii.gz'
	
	trim.inputs.in_file = imgFile 
	trim.inputs.out_file = aslVolFile
	trim.inputs.end_index = 1
	aslVol = trim.run()
	
	trim.inputs.out_file = pdVolFile
	trim.inputs.end_index = 2
	trim.inputs.begin_index = 1
	pdVol = trim.run()
	
	return aslVolFile
	return pdVolFile	

# Slice-timing correction


# NU correction

# M0 Calculation
# def M0calc(pdVol):
	
# Return corrected volume

# Run the defined functions
getImgMeta(img)
splitASLvols(imgFile, path)

