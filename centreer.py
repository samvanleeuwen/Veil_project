#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:02:42 2024

@author: samvanleeuwen
"""

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.io import fits  # to read in FITS files
import os  # os.path to manipulate file paths
import glob  # finding pathnames (to search for certain fits files in folders and subfolders)
import numpy as np  # math applied to arrays (important, no need to read pixel for pixel!)
from matplotlib import pyplot as plt  # plot library
from astropy.visualization import ZScaleInterval  # create minimum and maximum Z values for plotting
from photutils.datasets import load_star_image
from astropy.stats import sigma_clipped_stats
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from twirl import find_peaks
from photutils.aperture import CircularAperture
from astropy.wcs.utils import proj_plane_pixel_scales
from twirl import gaia_radecs
from twirl.geometry import sparsify
from twirl import compute_wcs
from reproject import reproject_interp
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import scipy.stats as stats

# Path to the FITS file to process
data_im = '/Users/samvanleeuwen/Downloads/Project Veil/Stack Ha alle datasets.fit'

try:
    # Open the FITS file and extract data
    with fits.open(data_im) as hdu_list:
        hdu = hdu_list[0]
        header = hdu.header
        data = hdu.data
        
    true_wcs = WCS(header)

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

except FileNotFoundError:
    print(f"File {data_im} not found.")
except ValueError as ve:
    print(f"ValueError: {ve}")
except Exception as e:
    print(f"An error occurred: {e}")

# Import list of fits files:
# Set path to calibration folder (path_in) and path to save files (path_out)
path_in = "/Users/samvanleeuwen/Downloads/Project Veil/DATA_skyview/"
path_out = "/Users/samvanleeuwen/Downloads/Project Veil/ouput_DATA/"

file_list = glob.glob(os.path.join(path_in, '*', '.fits*'), recursive=True)  # Search subfolders if recursive=True

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