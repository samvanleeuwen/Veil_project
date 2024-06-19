# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:58:45 2024

@author: milan
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
from skimage.draw import line
from matplotlib.patches import Polygon

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

# Background box coordinates
background_box = (2703, 3459, 2995, 3702)

def data_reduction(data, background_box):
    x1, y1, x2, y2 = background_box
    background_data = data[y1:y2, x1:x2]
    background_median = np.median(background_data)
    return data - background_median

def create_diagonal_box(x1, y1, box_width, box_height, angle):
    """
    Create a single diagonal box.
    
    Parameters:
    - x1, y1: Top-left coordinates of the box.
    - box_width: Width of the box.
    - box_height: Height of the box.
    - angle: Rotation angle of the box (in degrees).
    
    Returns:
    - x_coords, y_coords: Arrays of the x and y coordinates of the box corners.
    """
    angle_rad = np.deg2rad(angle)
    
    # Calculate the coordinates of the box corners
    x2 = x1 + box_width * np.cos(angle_rad)
    y2 = y1 + box_width * np.sin(angle_rad)
    x3 = x2 + box_height * np.cos(angle_rad + np.pi/2)
    y3 = y2 + box_height * np.sin(angle_rad + np.pi/2)
    x4 = x1 + box_height * np.cos(angle_rad + np.pi/2)
    y4 = y1 + box_height * np.sin(angle_rad + np.pi/2)
    
    # Create arrays of the coordinates
    x_coords = np.array([x1, x2, x3, x4, x1])
    y_coords = np.array([y1, y2, y3, y4, y1])
    
    return x_coords, y_coords

def calculate_intensity_profile(data, x_coords, y_coords, box_length):
    """
    Calculate the intensity profile along a diagonal box.
    
    Parameters:
    - data: 2D array of the image data.
    - x_coords, y_coords: Coordinates of the box corners.
    - box_length: Length of the diagonal box.
    
    Returns:
    - intensity_profile: Array of median intensity values along the diagonal box.
    """
    intensity_profile = []
    for i in range(box_length):
        # Calculate the coordinates of the line segment
        x_start = int(x_coords[0] + i * (x_coords[1] - x_coords[0]) / box_length)
        y_start = int(y_coords[0] + i * (y_coords[1] - y_coords[0]) / box_length)
        x_end = int(x_coords[3] + i * (x_coords[2] - x_coords[3]) / box_length)
        y_end = int(y_coords[3] + i * (y_coords[2] - y_coords[3]) / box_length)
        
        # Get the pixel values along the line segment
        rr, cc = line(y_start, x_start, y_end, x_end)
        line_data = data[rr, cc]
        
        # Calculate the median intensity value along the line segment
        median_intensity = np.median(line_data)
        intensity_profile.append(median_intensity)
    
    return np.array(intensity_profile)

def plot_intensity_profiles(data_list, boxes, background_box, box_length):
    """
    Plots the intensity profiles of diagonal boxes at specified positions.
    """
    colors = ['b', 'g', 'r', 'c', 'm']
    for box_idx, box in enumerate(boxes):
        plt.figure()
        for data_idx, (data, title) in enumerate(zip(data_list, titles)):
            reprojected_data, tgt_wcs = data
            reprojected_data = data_reduction(reprojected_data, background_box)
            x_coords, y_coords = box
            intensity_profile = calculate_intensity_profile(reprojected_data, x_coords, y_coords, box_length)
            normalized_intensity_profile = (intensity_profile - np.min(intensity_profile)) / (np.max(intensity_profile) - np.min(intensity_profile))
            plt.plot(np.arange(len(normalized_intensity_profile)), normalized_intensity_profile, color=colors[data_idx], label=f'{title}')
        
        plt.xlabel('Pixel index along diagonal')
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

# Example usage of create_diagonal_box
box_width = 500
box_height = 10
angle = 30  # Rotation angle of the box in degrees
box_length = 500  # Length of the diagonal box for intensity profile calculation

# Define top-left corner coordinates of the boxes
top_left_coords = [
    (1120, 2584),
    (2264, 304),
    (1608, 1648),
    (121, 2790)
]

# Generate diagonal boxes
diagonal_boxes = [create_diagonal_box(x, y, box_width, box_height, angle) for x, y in top_left_coords]

# Plot intensity profiles for each diagonal box
plot_intensity_profiles(reprojected_data_list, diagonal_boxes, background_box, box_length)

# Plot all reprojected FITS files with diagonal boxes
plt.figure(figsize=(18, 8))

for i, (data, title) in enumerate(zip(reprojected_data_list, titles)):
    reprojected_data, tgt_wcs = data
    ax = plt.subplot(1, len(fits_files), i + 1, projection=tgt_wcs)
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(reprojected_data)
    ax.imshow(reprojected_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    
    # Plot the diagonal boxes
    for x_coords, y_coords in diagonal_boxes:
        poly = Polygon(np.column_stack([x_coords, y_coords]), closed=True, edgecolor='r', facecolor='none', linewidth=1)
        ax.add_patch(poly)
    
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