# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:06:42 2024

@author: milan
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import os  # os.path to manipulate file paths
import glob  # finding pathnames (to search for certain fits files in folders and subfolders)



# Path to the FITS file to process
data_im = 'C:/Uva/Jaar 1 dingen/Project/Stacks project/Stack van alle datasets/GedownloadHa.fits'

try:
    # Open the FITS file and extract data
    with fits.open(data_im) as hdu_list:
        hdu = hdu_list[0]
        header = hdu.header
        data = hdu.data
        
    true_wcs = WCS(header)
    print(header)

    # Plot the FITS image with WCS coordinates
    plt.figure()
    #ax = fig.add_subplot(projection=true_wcs)

    # Apply ZScaleInterval to set the display range
    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)

    x_center = header['CRVAL1']
    y_center = header['CRVAL2']
    
    # Define the provided pixel scale in degrees
    pixel_scale = 0.431 / 3600

    left = x_center - (header['CRPIX1'] * pixel_scale)
    lower = y_center - (header['CRPIX2'] * pixel_scale)
    
    right = left + (4096*pixel_scale)
    top = lower + (4096*pixel_scale)
    
    extent_image = [left, right, lower, top]

    plt.imshow(data, vmin=z1, vmax=z2, cmap='viridis', origin='lower', extent=(left, right, lower, top))
    #plt.set_title('FITS Imagewith WCS Coordinates')
    
    plt.xlabel('RA')
    plt.ylabel('Dec')
    
    plt.ylim(top ,lower)
    plt.xlim(left, right)
    # Show plot
    plt.show()
    print(left)
    print(right)
    print(lower)
    print(top)

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

for getal in range(0,3):
    getal = getal + 1
    hdu_list= fits.open(file_list[getal])[0]
    header = hdu_list.header
    data, true_wcs = hdu_list.data, WCS(hdu_list.header)

    # Plot the FITS image with WCS coordinates
    plt.figure()
    #ax = fig.add_subplot(projection=true_wcs)

    # Apply ZScaleInterval to set the display range
    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)
    
    x_center = header['CRVAL1']
    y_center = header['CRVAL2']
    
    # Define the provided pixel scale in degrees
    pixel_scale = 1.76 / 3600

    left = x_center - (header['CRPIX1'] * pixel_scale)
    lower = y_center - (header['CRPIX2'] * pixel_scale)
    
    right = left + (1024*pixel_scale)
    top = lower + (1024*pixel_scale)
    
    extent_image = [left, right, lower, top]

    plt.imshow(data, vmin=z1, vmax=z2, cmap='viridis', origin='lower', extent=(left, right, lower, top))
    #plt.set_title('FITS Imagewith WCS Coordinates')
    
    plt.xlabel('RA')
    plt.ylabel('Dec')
    
    plt.ylim(top, lower)
    plt.xlim(right, left)

    plt.show()