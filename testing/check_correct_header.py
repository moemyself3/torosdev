import os

import pandas as pd
import numpy as np
import astropy.wcs as pywcs
import scipy.ndimage
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import twirl
from photutils.detection import DAOStarFinder
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats

from config import Configuration
from libraries.utils import Utils

path = Configuration.CLEAN_DIRECTORY + "2024-11-02/FIELD_2b.022/"
files = os.listdir(path)
files.sort()
file = files[0]
maskpath = Configuration.CALIBRATION_DIRECTORY + "mask.fits"
mask = fits.getdata(maskpath).astype('bool')

def correct_header(filename=path+file, mask=mask):
    img, header = fits.getdata(filename , header=True)

    center = SkyCoord(Configuration.RA, Configuration.DEC, unit=['deg','deg'])
    pixel = Configuration.PIXEL_SIZE * u.arcsec
    fov = np.max(np.shape(img)) * pixel.to(u.deg)

    all_stars = twirl.gaia_radecs(center, 1.25 * fov)
    all_stars_sparse = twirl.geometry.sparsify(all_stars, 0.01)[0:30]

    mean, median, std = sigma_clipped_stats(img, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=100, exclude_border=True, peakmax=45000)
    sources = daofind(img - median, mask=mask)

    return sources, all_stars, all_stars_sparse

