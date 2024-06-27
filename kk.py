

import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

# Open some FITS image
hdu = fits.open("...")[0]

# get the center of the image
ra, dec = hdu.header["RA"], hdu.header["DEC"]
center = SkyCoord(ra, dec, unit=["deg", "deg"])

# and the size of its field of view
pixel = 0.43 * u.arcsec  # known pixel scale
shape = hdu.data.shape
fov = np.max(shape) * pixel.to(u.deg)