
# coding: utf-8

# In[3]:

get_ipython().magic(u'matplotlib')
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib


# In[13]:

def plot_slice(fname, offset, vmin, vmax):
    # Load the image and collect the data orientation information
    img = nib.load(fname)
    data = img.get_data()
    aff = img.affine
    # Find the center of the brain matrix
    ctr = np.dot(np.linalg.inv(aff), [0, 0, 0, 1])[:3]
    # Plot the data
    plt.imshow(np.rot90(data[:, :,ctr[2] + offset]), cmap="viridis", vmin=vmin, vmax=vmax)
    plt.gca().set_axis_off()


# In[14]:

plot_slice("/Users/catcamacho/Documents/asl_practice/ko/wrcbf_gm.nii.gz",1,0,200)

