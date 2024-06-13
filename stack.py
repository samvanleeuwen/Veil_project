# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 17:40:00 2024

@author: milan
"""

from astropy.wcs import WCS
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import os  # os.path to manipulate file paths
import glob  # finding pathnames (to search for certain fits files in folders and subfolders)

# Path to the FITS file to process
data_im = 'C:/Uva/Jaar 1 dingen/Project/Stacks project/Stack van alle datasets/GedownloadO3.fits'

try:
    # Open the FITS file and extract data
    with fits.open(data_im) as hdu_list:
        hdu = hdu_list[0]
        header = hdu.header
        data = hdu.data
 
    wcs = WCS(header)

    wx1, wy1 = wcs.wcs_world2pix(313.9,30.87, 1)
    wx2, wy2 = wcs.wcs_world2pix(314.5,31.34, 1)
    
    
    fig = plt.figure()
    ax = plt.subplot(111, projection=wcs)
    
    # Apply ZScaleInterval to set the display range
    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)
    
    ax.imshow(data, cmap='gray', vmin=z1, vmax=z2, interpolation=None, origin='lower')
    ax.set_xlabel("Right Ascension [degrees]")
    ax.set_ylabel("Declination [degrees]")
    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    ax.set_autoscale_on(False)
    ax.set_xlim(wx1,wx2)
    ax.set_ylim(wy1,wy2)
    
    plt.plot()
    
except FileNotFoundError:
    print(f"File {data_im} not found.")
except ValueError as ve:
    print(f"ValueError: {ve}")
except Exception as e:
    print(f"An error occurred: {e}")
    
# Import list of fits files:
# Set path to calibration folder (path_in) and path to save files (path_out)
path_in = "C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen"
path_out = "C:/Uva/Jaar 1 dingen/Project/Stacks project/Stacks gecalibreerd met x-ray en röntgen/"

file_list = glob.glob(os.path.join(path_in, '**', '*.fits*'), recursive=True)  # Search subfolders if recursive=True

# Debugging output to check file_list contents
print(f"Files found: {file_list}")

sorted_list = sorted(file_list, key=os.path.getmtime)

# Debugging output to check sorted_list contents
print(f"Sorted files: {sorted_list}")

for getal in range(0,5):
    getal = getal + 1
    hdu_list= fits.open(file_list[getal])[0]
    header = hdu_list.header
    data, true_wcs = hdu_list.data, WCS(hdu_list.header)

    wcs = WCS(header)

    wx1, wy1 = wcs.wcs_world2pix(313.9,30.87, 1)
    wx2, wy2 = wcs.wcs_world2pix(314.5,31.34, 1)
    
    
    fig = plt.figure()
    ax = plt.subplot(111, projection=wcs)
    
    # Apply ZScaleInterval to set the display range
    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)
    
    ax.imshow(data, cmap='gray', vmin=z1, vmax=z2, interpolation=None, origin='lower')
    ax.set_xlabel("Right Ascension [degrees]")
    ax.set_ylabel("Declination [degrees]")
    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    ax.set_autoscale_on(False)
    ax.set_xlim(wx1,wx2)
    ax.set_ylim(wy1,wy2)
    
    plt.plot()