import os
import numpy as np
import nipy as ni
import nipype as npe
import nibabel as nib
from glob import glob
from nipype.interfaces.nipy.preprocess import Trim
from nipype.interfaces import fsl

trim = Trim()
st = fsl.SliceTimer()

path = '/Users/catcamacho/Documents/asl_practice/210-NART2'
img_file = path + '/raw/asl.nii.gz'
asl_vol_file = path + '/aslVol.nii.gz'
pd_vol_file = path + '/pdVol.nii.gz'
st_corr_vol_file = path + '/st_aslVol.nii.gz'
n3_corr_vol_file = path + '/nust_aslVol.nii.gz'
img_data = nib.load(img_file)


def data_type(img_file):
    """
    Determine file type and transform to NIFTI if need be.
    """
    pass


def get_img_meta(img):
    """
    Collect file properties to be used as sanity checks down the line.
    """
    header = img.header
    dim = header.get_data_shape()
    matrix = dim[0:2]
    numslices = dim[2]
    numvols = dim[3]
    dimensions = {'matrix': matrix, 'numslices': numslices, 'numvols': numvols}

    print dimensions
    return dimensions


def split_asl_vols(img_file, asl_vol_file, pd_vol_file):
    """
    Read in the combo nifti and and split the ASL and PD Volumes.
    """
    trim.inputs.in_file = img
    trim.inputs.out_file = asl_vol_file
    trim.inputs.end_index = 1
    trim.run()
    trim.inputs.out_file = pd_vol_file
    trim.inputs.end_index = 2
    trim.inputs.begin_index = 1
    trim.run()


def slicetime_asl(asl_vol_file, st_corr_vol_file):
    """
    Correct slice timing using FSL-need to replace with a custom one.
    """
    st.inputs.time_repetition = 4.674
    st.inputs.slice_direction = 3
    st.inputs.interleaved = False
    st.inputs.in_file = asl_vol_file
    st.inputs.out_file = st_corr_vol
    st.inputs.output_type = 'NIFTI_GZ'
    st.run()


def nonuniformity_correction_asl(st_corr_vol_file):
    """
    Correct nonuniformity slice by slice low pass filtering via 3D FFT
    """
    pass


def downsample_anat():
    """
    Downsample anat to 2mm cubic voxels (FS/iBEAT outputs so no corrections)
    Also downsample aparc+aseg parcellation (if applicable; aseg and desikan
    killany atlas), the WM seg, and the GM seg.
    """
    pass


def upsample_asl():
    """
    Upsample the ASL volume interpolating the neighboring values. Final
    voxel size will be 2mm cubed.
    """
    pass


def coregister_vols():
    """
    Coregister the two inputed volumes. Apply same transformation matrix
    to the segmentation volumes.
    """
    pass


def calculate_m0(pd_vol):
    """
    Calculate M0 value from the PD volume.
    """
    pd = nib.load(pd_vol)


def main():
    get_img_meta(img)
    split_asl_vols(img_file, asl_vol_file, pd_vol_file)
    slicetime_asl(asl_vol_file, st_corr_vol_file)
    nonuniformity_correction_asl(st_corr_vol_file)
    calculate_m0(pd_vol)

if __name__ == '__main__':
    main()
