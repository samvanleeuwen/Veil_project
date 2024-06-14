# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:58:45 2024

@author: milan
"""

import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
from astropy.visualization import ZScaleInterval
from reproject import reproject_interp
from astropy.visualization import ZScaleInterval

filename = 'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/GedownloadHa.fits'

# Path to the FITS files
reference_fits = 'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/GedownloadHa.fits'  # FITS file with the desired WCS
target_fits = 'C:/Uva/Jaar 1 dingen/Project/Stacks project/Filters x-ray en röntgen/GedownloadWise12.fits'     # FITS file to be aligned

# Open the reference FITS file
with fits.open(reference_fits) as ref_hdu_list:
    ref_hdu = ref_hdu_list[0]
    ref_data = ref_hdu.data
    ref_wcs = WCS(ref_hdu.header)

# Open the target FITS file
with fits.open(target_fits) as tgt_hdu_list:
    tgt_hdu = tgt_hdu_list[0]
    tgt_data = tgt_hdu.data
    tgt_wcs = WCS(tgt_hdu.header)

from reproject import reproject_interp

# Reproject the target image to the WCS of the reference image
reprojected_data, footprint = reproject_interp((tgt_data, tgt_wcs), ref_wcs, shape_out=ref_data.shape)

# Apply ZScaleInterval to set the display range
z = ZScaleInterval()
ref_z1, ref_z2 = z.get_limits(ref_data)
tgt_z1, tgt_z2 = z.get_limits(reprojected_data)

# Plot the reference image
fig, ax = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ref_wcs})
ax[0].imshow(ref_data, vmin=ref_z1, vmax=ref_z2, cmap='gray', origin='lower')
ax[0].set_title('Reference Image')
ax[0].coords.grid(color='white', ls='solid')
ax[0].set_xlabel('RA')
ax[0].set_ylabel('Dec')

# Plot the reprojected target image
ax[1].imshow(reprojected_data, vmin=tgt_z1, vmax=tgt_z2, cmap='gray', origin='lower')
ax[1].set_title('Reprojected Target Image')
ax[1].coords.grid(color='white', ls='solid')
ax[1].set_xlabel('RA')
ax[1].set_ylabel('Dec')

plt.show()