#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:39:57 2024

@author: heikamp


start up script for astronomy image analysis with python

"""
# Import packages
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



# Import fits file
data_im = 'C:/Uva/Jaar 1 dingen/Project/Stacks project/Stack van alle datasets/Stack Ha alle datasets.fit'

try:
    hdu = fits.open(data_im)[0]
    header = hdu.header
    data = hdu.data
    data, true_wcs = hdu.data, WCS(hdu.header)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    xy = find_peaks(data)[0:20]
    
    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)
    plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')  # Display the main image using ZScaleInterval
    plt.imshow(data, vmin=np.median(data), vmax=3 * np.median(data), cmap="gray", alpha=0.5)  # Overlay the background image
    _ = CircularAperture(xy, r=10.0).plot(color="y")
    fov = (data.shape * proj_plane_pixel_scales(true_wcs))[0]
    center = true_wcs.pixel_to_world(*np.array(data.shape) / 2)
    all_radecs = gaia_radecs(center, 1.2 * fov)

    # we only keep stars 0.01 degree apart from each other
    all_radecs = sparsify(all_radecs, 0.01)
    # we only keep the 12 brightest stars from gaia
    wcs = compute_wcs(xy, all_radecs[0:30], tolerance=10)
    # plotting to check the WCS
    radecs_xy = np.array(wcs.world_to_pixel_values(all_radecs))
    plt.imshow(data, vmin=np.median(data), vmax=3 * np.median(data), cmap="gray", alpha=0.5)
    _ = CircularAperture(radecs_xy, 5).plot(color="y", alpha=0.5)
    plt.show()
except FileNotFoundError:
    print(f"File {data_im} not found.")
except Exception as e:
    print(f"An error occurred while opening the FITS file: {e}")

# Import list of fits files:
# Set path to calibration folder (path_in) and path to save files (path_out)
path_in = "C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen"
path_out = "C:/Uva/Jaar 1 dingen/Project/Stacks project/Stacks gecalibreerd met x-ray en röntgen/"

file_list = glob.glob(os.path.join(path_in, '**', '*.fit*'), recursive=True)  # Search subfolders if recursive=True

# Debugging output to check file_list contents
print(f"Files found: {file_list}")

sorted_list = sorted(file_list, key=os.path.getmtime)

# Debugging output to check sorted_list contents
print(f"Sorted files: {sorted_list}")

for getal in range(0,3):
    hdu_list= fits.open(sorted_list[getal])[0]
    header = hdu_list.header
    data, true_wcs = hdu_list.data, WCS(hdu_list.header)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    xy = find_peaks(data)[0:20]

    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)
    plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')
    plt.imshow(data, vmin=np.median(data), vmax=3 * np.median(data), cmap="grey")
    _ = CircularAperture(xy, r=10.0).plot(color="y")
    fov = (data.shape * proj_plane_pixel_scales(true_wcs))[0]
    center = true_wcs.pixel_to_world(*np.array(data.shape) / 2)
    all_radecs = gaia_radecs(center, 1.2 * fov)

    # we only keep stars 0.01 degree apart from each other
    all_radecs = sparsify(all_radecs, 0.01)
    # we only keep the 12 brightest stars from gaia
    wcs = compute_wcs(xy, all_radecs[0:30], tolerance=10)
    # plotting to check the WCS
    radecs_xy = np.array(wcs.world_to_pixel_values(all_radecs))
    plt.imshow(data, vmin=np.median(data), vmax=3 * np.median(data), cmap="grey", alpha=0.5)
    _ = CircularAperture(radecs_xy, 5).plot(color="y", alpha=0.5)
    plt.show()

if len(sorted_list) > 1:
    try:
        hdu_list = fits.open(sorted_list[2])[0]  # 1 = 2nd item in list (index starts at 0)
        header = hdu_list.header
        data = hdu_list.data

        # Simple manipulation and save fits file
        data_transposed = data.T
        data_sub = data_transposed - data_transposed  # Example manipulation
        data_median = np.median(data)

        hdu = fits.PrimaryHDU(data_sub, header=header)
        hdu.writeto(os.path.join(path_out, "filename.fits"), overwrite=True)
    except Exception as e:
        print(f"An error occurred while processing the FITS file: {e}")
else:
    print("No FITS files found or not enough files in the specified directory.")