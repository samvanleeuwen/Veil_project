# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:58:45 2024

@author: milan
"""

import numpy as np
from astropy.io import fits
from reproject import reproject_interp
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
from skimage.draw import line
from matplotlib.patches import Polygon

from astropy.wcs import WCS

plt.rcParams['text.usetex'] = False





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
    angle_rad = np.deg2rad(angle)
    x2 = x1 + box_width * np.cos(angle_rad)
    y2 = y1 + box_width * np.sin(angle_rad)
    x3 = x2 + box_height * np.cos(angle_rad + np.pi/2)
    y3 = y2 + box_height * np.sin(angle_rad + np.pi/2)
    x4 = x1 + box_height * np.cos(angle_rad + np.pi/2)
    y4 = y1 + box_height * np.sin(angle_rad + np.pi/2)
    x_coords = np.array([x1, x2, x3, x4, x1])
    y_coords = np.array([y1, y2, y3, y4, y1])
    return x_coords, y_coords

def calculate_intensity_profile(data, x_coords, y_coords, box_length):
    intensity_profile = []
    for i in range(box_length):
        x_start = int(x_coords[0] + i * (x_coords[3] - x_coords[0]) / box_length)
        y_start = int(y_coords[0] + i * (y_coords[3] - y_coords[0]) / box_length)
        x_end = int(x_coords[1] + i * (x_coords[2] - x_coords[1]) / box_length)
        y_end = int(y_coords[1] + i * (y_coords[2] - y_coords[1]) / box_length)
        rr, cc = line(y_start, x_start, y_end, x_end)
        line_data = data[rr, cc]
        median_intensity = np.median(line_data)
        intensity_profile.append(median_intensity)
    return np.array(intensity_profile)

def plot_intensity_profiles(data_list, boxes, background_box, box_length):
    colors = ['r', 'g', 'b', 'c', 'm']
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

with fits.open(fits_files[0]) as ref_hdu_list:
    ref_hdu = ref_hdu_list[0]
    ref_data = ref_hdu.data.astype(float)

reprojected_data_list = []

for file, title in zip(fits_files, titles):
    with fits.open(file) as tgt_hdu_list:
        tgt_hdu = tgt_hdu_list[0]
        tgt_data = tgt_hdu.data.astype(float)
        reprojected_data, footprint = reproject_interp((tgt_data, WCS(tgt_hdu.header)), WCS(ref_hdu.header), shape_out=ref_data.shape)
        reprojected_data_list.append((reprojected_data, WCS(tgt_hdu.header)))

box_width = 14
box_height = 150
angle = 20.1
box_length = 150

top_left_coords = [
    (2670, 2425),
    (1395, 2660)
]

diagonal_boxes = [create_diagonal_box(x, y, box_width, box_height, angle) for x, y in top_left_coords]

plot_intensity_profiles(reprojected_data_list, diagonal_boxes, background_box, box_length)

fig, axes = plt.subplots(1, len(fits_files), figsize=(18, 8))

def degrees_to_hms(degrees):
    total_seconds = degrees * 3600
    hours = int(total_seconds // 3600)
    total_seconds %= 3600
    minutes = int(total_seconds // 60)
    seconds = total_seconds % 60
    return f"{hours}h {minutes}m"

def degrees_to_dms(degrees):
    total_seconds = degrees * 3600
    degrees = int(total_seconds // 3600)
    total_seconds %= 3600
    minutes = int(total_seconds // 60)
    seconds = total_seconds % 60
    return f"{degrees}° {minutes}\'"

getal = 0

for i, (data, title) in enumerate(zip(reprojected_data_list, titles)):
    getal = getal + 1
    ax = axes[i]
    reprojected_data, tgt_wcs = data
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(reprojected_data)
    ax.imshow(reprojected_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    
    if getal==2:
        
        # Plot the diagonal boxes
        for x_coords, y_coords in diagonal_boxes:
            poly = Polygon(np.column_stack([x_coords, y_coords]), closed=True, edgecolor='g', facecolor='none', linewidth=1)
            ax.add_patch(poly)
    
        
        # Add background box
        x1, y1, x2, y2 = background_box
        if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
            rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=1, edgecolor='b', facecolor='none')
            ax.add_patch(rect)
        else:
            print(f"Warning: Background box {background_box} exceeds image dimensions for {title}")
    if getal == 1:
                
        # Plot the diagonal boxes
        for x_coords, y_coords in diagonal_boxes:
            poly = Polygon(np.column_stack([x_coords, y_coords]), closed=True, edgecolor='r', facecolor='none', linewidth=1)
            ax.add_patch(poly)
    
        
        # Add background box
        x1, y1, x2, y2 = background_box
        if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
            rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=1, edgecolor='b', facecolor='none')
            ax.add_patch(rect)
        else:
            print(f"Warning: Background box {background_box} exceeds image dimensions for {title}")

    if i == 0 or i == 1:
        color = 'y' if i == 0 else 'y'
        for box_idx, (x_coords, y_coords) in enumerate(diagonal_boxes):
            poly = Polygon(np.column_stack([x_coords, y_coords]), closed=True, edgecolor=color, facecolor='none', linewidth=1, alpha=0.5)
            ax.add_patch(poly)
            # Adding labels for the boxes outside the boxes
            mid_x = (x_coords[0] + x_coords[2]) / 2
            mid_y = (y_coords[0] + y_coords[2]) / 2
            offset = 40  # Offset to move the text outside the box
            ax.text(mid_x + offset, mid_y + offset, str(box_idx + 1), color=color, fontsize=11, ha='center', va='center', alpha=0.7)
        x1, y1, x2, y2 = background_box
        if x2 <= reprojected_data.shape[1] and y2 <= reprojected_data.shape[0]:
            rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, linewidth=1, edgecolor='b', facecolor='none')
            ax.add_patch(rect)
        else:
            print(f"Warning: Background box {background_box} exceeds image dimensions for {title}")

    ax.set_title(f'{title}')
    ax.set_xlabel('RA')
    ax.set_xlim(0, reprojected_data.shape[1])
    ra_ticks = np.linspace(314.04326958715717, 314.53365180937936, 5)
    ax.set_xticks(np.linspace(0, reprojected_data.shape[1], 5))
    ax.set_xticklabels([degrees_to_hms(tick) for tick in ra_ticks], rotation=45, ha='right')

    if i == 0:
        ax.set_ylabel('Dec')
        ax.set_ylim(0, reprojected_data.shape[0])
        dec_ticks = np.linspace(30.85197440095892, 31.342356623181143, 5)
        ax.set_yticks(np.linspace(0, reprojected_data.shape[0], 5))
        ax.set_yticklabels([degrees_to_dms(tick) for tick in dec_ticks])
    else:
        ax.set_ylabel('')
        ax.yaxis.set_tick_params(labelleft=False)

plt.subplots_adjust(wspace=0.4, hspace=0.1)
plt.show()