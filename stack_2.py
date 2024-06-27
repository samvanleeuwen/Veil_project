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

try:
    hdu = fits.open(data_im)[0]
    header = hdu.header
    data = hdu.data

    print(header)

    z = ZScaleInterval()
    z1, z2 = z.get_limits(data)
    plt.imshow(data, vmin=z1, vmax=z2, cmap='gray')
    plt.show()
except FileNotFoundError:
    print(f"File {data_im} not found.")
except Exception as e:
    print(f"An error occurred while opening the FITS file: {e}")

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

file_list = glob.glob(os.path.join(path_in, '**', '*.fit*'), recursive=True)  # Search subfolders if recursive=True

# Debugging output to check file_list contents
print(f"Files found: {file_list}")

sorted_list = sorted(file_list, key=os.path.getmtime)

# Debugging output to check sorted_list contents
print(f"Sorted files: {sorted_list}")

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
