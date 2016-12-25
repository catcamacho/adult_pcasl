import os
import numpy as np
import nibabel as nib
from glob import glob
import nipy as ni
import nipype as npe
from nipype.interfaces.nipy.preprocess import Trim
trim = Trim()
from nipype.interfaces import fsl
st = fsl.SliceTimer()

path = '/Users/catcamacho/Documents/asl_practice/210-NART2'
imgFile = path + '/raw/asl.nii.gz'
img = nib.load(imgFile)
aslVolFile = path + '/aslVol.nii.gz'
pdVolFile = path + '/pdVol.nii.gz'
	

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
def splitASLvols(imgFile, aslVolFile, pdVolFile):	
	trim.inputs.in_file = imgFile 
	trim.inputs.out_file = aslVolFile
	trim.inputs.end_index = 1
	trim.run()
	
	trim.inputs.out_file = pdVolFile
	trim.inputs.end_index = 2
	trim.inputs.begin_index = 1
	trim.run()

# Slice-timing correction
def sliceTimeASL(aslVolFile, path):
	st.inputs.time_repetition = 4.674
	st.inputs.slice_direction = 3
	st.inputs.interleaved = False
	st.inputs.in_file = aslVolFile
	st.inputs.out_file = path + '/st_aslVol.nii.gz'
	st.inputs.output_type = 'NIFTI_GZ'
	st.run()
	
# NU correction


# M0 Calculation
# def M0calc(pdVol):
	
# Return corrected volume

# Run the defined functions
getImgMeta(img)
splitASLvols(imgFile, aslVolFile, pdVolFile)
sliceTimeASL(aslVolFile, path)

