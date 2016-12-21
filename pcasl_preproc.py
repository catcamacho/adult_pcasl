import os
import numpy as np
import nibabel as nib
import nipy as nip

np.set_printoptions(precision=4, suppress=True)
path = '/Users/catcamacho/Box Sync/CBFBIRN_practice/210-NART2'
img = nib.load((path + '/raw/asl.nii.gz'))

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

#	return perf pd
#def splitVols(img,numVols):
#	aslVol = nib.nifti1.Nifti1Image(img.dataobj[:1]) #helpppp
#	pdVol = nib.nifti1.Nifti1Image(img.data[1:], img.header) #helpppp
#	aslVol.save((path + '/aslVol.nii.gz')) #helpppp
#	pdVol.save((path + '/pdVol.nii.gz')) #helpppp
#	return pdVol
#	return aslVol	

# Slice-timing correction

# NU correction

# M0 Calculation

# Return corrected volume

# Run the defined functions
getImgMeta(img)

