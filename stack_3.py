#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:46:48 2024

@author: samvanleeuwen
"""

import numpy as np
import math
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt

# Paths to the FITS files
fits_files = [
    '/Users/samvanleeuwen/Downloads/Project Veil/Stack Ha alle datasets.fit',
    '/Users/samvanleeuwen/Downloads/Project Veil/Stack OIII alle datasets.fit',
]

# Titles for each image
titles = [
    'Ha',
    'OIII'
]

# Define a function to create a single diagonal box
def create_diagonal_box(image_width, image_height, box_width, box_height, angle):
    # Convert angle to radians
    angle_rad = np.deg2rad(angle)
    
    # Calculate the coordinates of the box corners
    x1, y1 = 100, 100
    x2 = box_width * np.cos(angle_rad)
    y2 = box_width * np.sin(angle_rad)
    x3 = x2 + box_height * np.cos(angle_rad + np.pi / 2)
    y3 = y2 + box_height * np.sin(angle_rad + np.pi / 2)
    x4 = x1 + box_height * np.cos(angle_rad + np.pi / 2)
    y4 = y1 + box_height * np.sin(angle_rad + np.pi / 2)
    
    # Ensure the box coordinates are within the image dimensions
    x = np.array([x1, x2, x3, x4, x1])
    y = np.array([y1, y2, y3, y4, y1])
    
    # Translate box to the center of the image
    x += image_width // 2
    y += image_height // 2
    
    return [(x, y)]

# Example usage
image_width = 4096  # Positive width
image_height = 4096  # Positive height
box_width = 1000
box_height = 100
angle = 240  # Angle of the diagonal box in degrees

boxes = create_diagonal_box(image_width, image_height, box_width, box_height, angle)

# Define a function to plot the intensity profiles for each box
def plot_intensity_profiles(data_list, boxes):
    """
    Plots the intensity profiles of vertical boxes at specified positions.
    """
    colors = ['b', 'g', 'r', 'c', 'm']
    for box_idx, (x, y) in enumerate(boxes):
        plt.figure(figsize=(12, 6))
        
        for data_idx, (data, title) in enumerate(zip(data_list, titles)):
            reprojected_data, tgt_wcs = data
            
            # Ensure the box coordinates are within the bounds of the image
            if np.any(x < 0) or np.any(y < 0) or np.any(x > reprojected_data.shape[1]) or np.any(y > reprojected_data.shape[0]):
                print(f"Warning: Box {box_idx + 1} exceeds image dimensions for {title}")
                continue
            
            # Create a mask for the box
            mask = np.zeros(reprojected_data.shape, dtype=bool)
            rr, cc = np.array(np.round(y), dtype=int), np.array(np.round(x), dtype=int)
            mask[rr, cc] = True
            
            # Extract the intensity profile
            box_data = reprojected_data[mask]
            
            if box_data.size == 0:
                print(f"Warning: Box {box_idx + 1} results in an empty slice for {title}")
                continue
            
            # Calculate the intensity profile along the x-axis (horizontal direction)
            intensity_profile = np.mean(box_data, axis=0)
            
            # Normalize the intensity profile
            if np.ptp(intensity_profile) == 0:  # Check if all values are the same
                normalized_intensity_profile = intensity_profile
            else:
                normalized_intensity_profile = (intensity_profile - np.min(intensity_profile)) / (np.max(intensity_profile) - np.min(intensity_profile))
            
            # Plotting
            x_positions = np.arange(x[0], x[-1] + 1)  # Ensure x_positions cover the correct range
            plt.plot(x_positions, normalized_intensity_profile, color=colors[data_idx], label=title)
        
        plt.xlabel('Pixel x-coordinate')
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
    
    for box_idx, (x, y) in enumerate(boxes):
        if np.all(x >= 0) and np.all(y >= 0) and np.all(x <= reprojected_data.shape[1]) and np.all(y <= reprojected_data.shape[0]):
            polygon = plt.Polygon(np.column_stack([x, y]), linewidth=1, edgecolor='r', facecolor='none')
            ax.add_patch(polygon)
        else:
            print(f"Warning: Box {box_idx + 1} exceeds image dimensions for {title}")
    
    ax.set_title(f'{title} with Analyzed Boxes')
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    
plt.tight_layout()
plt.show()

