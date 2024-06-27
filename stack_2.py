#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 15:02:42 2024

@author: samvanleeuwen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:39:57 2024

@author: heikamp

Start-up script for astronomy image analysis with Python
"""

# Import packages
<<<<<<< HEAD
from astropy.io import fits  # to read in FITS files
import os  # os.path to manipulate file paths
import glob  # finding pathnames (to search for certain fits files in folders and subfolders)
import numpy as np  # math applied to arrays (important, no need to read pixel for pixel!)
from matplotlib import pyplot as plt  # plot library
from astropy.visualization import ZScaleInterval  # create minimum and maximum Z values for plotting
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import twirl

# Import fits file
data_im = 'C:/Users/elwin/.conda/20240518/H-a stacked'
=======
# Paths to the FITS files
fits_files = [
    'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/GedownloadHa.fits',
    'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/GedownloadO3.fits',
    'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/DSS1 2.fits',
    'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/Wise 12 2.fits',
    'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/PSPC2.0 2.fits'
]

# Titles for each image
titles = [
    'Ha',
    'OIII',
    'DSS 1 red',
    'Wise 12',
    'PSPC 2.0'
]
>>>>>>> 8e384777cc26612d6ad54ca4dccc731b6ff555ff

# Fixed Box coordinates [(x1, y1, x2, y2)]
boxes = [
    (1120, 2584, 1936, 2912),
    (2264, 304, 2696, 576),
    (1608, 1648, 2040, 1944),
    (121, 2790, 780, 3321)
]

background_box = (2703, 3459, 2995, 3702)

def data_reduction(data, background_box):
    x1, y1, x2, y2 = background_box
    background_data = data[y1:y2, x1:x2]
    background_median = np.median(background_data)
    return data - background_median

<<<<<<< HEAD
# Import list of fits files:
# Set path to calibration folder (path_in) and path to save files (path_out)
path_in = "C:/Users/elwin/.conda/20240518/Project_00"
path_out = "C:/Users/elwin/.conda/20240518/Project_01"

# Open some FITS image
hdu = fits.open(data_im)[0]

# get the center of the image
ra, dec = hdu.header["RA"], hdu.header["DEC"]
center = SkyCoord(ra, dec, unit=["deg", "deg"])

# and the size of its field of view
pixel = 0.66 * u.arcsec  # known pixel scale
shape = hdu.data.shape
fov = np.max(shape) * pixel.to(u.deg)
=======
# Define a function to plot the intensity profiles for each box
def plot_intensity_profiles(data_list, boxes, background_box):
    """
    Plots the intensity profiles of vertical boxes at specified positions.
    """
    colors = ['b', 'g', 'r', 'c', 'm']
    for box_idx, box in enumerate(boxes):
        plt.figure()
        for data_idx, (data, title) in enumerate(zip(data_list, titles)):
            reprojected_data, tgt_wcs = data
            print('oude')
            print(reprojected_data)
            reprojected_data = data_reduction(reprojected_data, background_box)
            print('nieuwe')
            print(reprojected_data)
            x1, y1, x2, y2 = box
            if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
                box_data = reprojected_data[y1:y2, x1:x2]
                intensity_profile = np.median(box_data, axis=1)  
                normalized_intensity_profile = (intensity_profile - np.min(intensity_profile)) / (np.max(intensity_profile) - np.min(intensity_profile))
                y_positions = np.arange(y1, y2)
                plt.plot(y_positions, normalized_intensity_profile, color=colors[data_idx], label=f'{title}')
            else:
                print(f"Warning: Box {box} exceeds image dimensions for {title}")
        
        plt.xlabel('Pixel y-coordinate')
        plt.ylabel('Normalized intensity')
        plt.title(f'Intensity Profiles for Box {box_idx + 1}')
        plt.legend()
        plt.show()
>>>>>>> 8e384777cc26612d6ad54ca4dccc731b6ff555ff

# Open the reference FITS file (first file in the list)
with fits.open(fits_files[0]) as ref_hdu_list:
    ref_hdu = ref_hdu_list[0]
    ref_data = ref_hdu.data.astype(float)  # Ensure data is float
    ref_wcs = WCS(ref_hdu.header)

# Initialize a list to store the reprojected data
reprojected_data_list = []

# Reproject each FITS file to the WCS of the reference image
for file, title in zip(fits_files, titles):
    with fits.open(file) as tgt_hdu_list:
        tgt_hdu = tgt_hdu_list[0]
        tgt_data = tgt_hdu.data.astype(float)  # Ensure data is float
        tgt_wcs = WCS(tgt_hdu.header)
        
        reprojected_data, footprint = reproject_interp((tgt_data, tgt_wcs), ref_wcs, shape_out=ref_data.shape)
        reprojected_data_list.append((reprojected_data, tgt_wcs))

# Plot intensity profiles for each box
plot_intensity_profiles(reprojected_data_list, boxes, background_box)

# Plot all reprojected FITS files with boxes
plt.figure(figsize=(18, 8))

for i, (data, title) in enumerate(zip(reprojected_data_list, titles)):
    reprojected_data, tgt_wcs = data
    ax = plt.subplot(1, len(fits_files), i + 1, projection=tgt_wcs)
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(reprojected_data)
    ax.imshow(reprojected_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    
    for box in boxes:
        x1, y1, x2, y2 = box
        if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
            rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=1, edgecolor='r', facecolor='none')
            ax.add_patch(rect)
        else:
            print(f"Warning: Box {box} exceeds image dimensions for {title}")

    # Add background box
    x1, y1, x2, y2 = background_box
    if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
        rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=1, edgecolor='g', facecolor='none')
        ax.add_patch(rect)
    else:
        print(f"Warning: Background box {background_box} exceeds image dimensions for {title}")
    
    ax.set_title(f'{title} with Analyzed Boxes')
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    
plt.tight_layout()
plt.show()
