# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:58:45 2024

@author: milan
"""

import numpy as np
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt

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

# Fixed Box coordinates [(x1, y1, x2, y2)]
boxes = [
    (1120, 2584, 1936, 2912),
    (2264, 304, 2696, 576),
    (1608, 1648, 2040, 1944),
    (121, 2790, 780, 3321)
]

# Define a function to plot the intensity profiles for each box
def plot_intensity_profiles(data_list, boxes):
    """
    Plots the intensity profiles of vertical boxes at specified positions.
    """
    colors = ['b', 'g', 'r', 'c', 'm']
    for box_idx, box in enumerate(boxes):
        plt.figure(figsize=(12, 6))
        
        for data_idx, (data, title) in enumerate(zip(data_list, titles)):
            reprojected_data, tgt_wcs = data
            x1, y1, x2, y2 = box
            if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
                box_data = reprojected_data[y1:y2, x1:x2]
                intensity_profile = np.mean(box_data, axis=1)
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
plot_intensity_profiles(reprojected_data_list, boxes)

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
    
    ax.set_title(f'{title} with Analyzed Boxes')
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    
plt.tight_layout()
plt.show()