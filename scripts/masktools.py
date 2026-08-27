import ccdproc
import astropy.units as u
from astropy.nddata import CCDData

import os

from libraries.preprocessing import Preprocessing

#path = "/Volumes/datadrive/coldharbor/flats/2024-11-02/"

#files = os.listdir(path)
#files.sort()

#flat1 = files[0]
#flat2 = files[-1]

filename = "/Volumes/datadrive/coldharbor/calibration/flat.fits"

ccd = CCDData.read(filename, unit="adu")
#ccd1 = CCDData.read(path+flat1, unit="adu")
#ccd2 = CCDData.read(path+flat2, unit="adu")

#ccd1.data, ccd1.header = Preprocessing.subtract_scaled_bias_dark(ccd1.data, ccd1.header)
#ccd2.data, ccd2.header = Preprocessing.subtract_scaled_bias_dark(ccd2.data, ccd2.header)

ccd.data, ccd.header = Preprocessing.clip_image(ccd.data, ccd.header)

#ccd1.data, ccd1.header = Preprocessing.clip_image(ccd1.data, ccd1.header)
#ccd2.data, ccd2.header = Preprocessing.clip_image(ccd2.data, ccd2.header)

ccd.header["NAXIS1"], ccd.header["NAXIS2"] = ccd.data.shape

#ccd1.header["NAXIS1"], ccd1.header["NAXIS2"] = ccd1.data.shape
#ccd2.header["NAXIS1"], ccd2.header["NAXIS2"] = ccd2.data.shape

#ratio = ccd2.divide(ccd1)

#print(f"Ratio mean: {ratio.data.mean()}")



maskr = ccdproc.ccdmask(ccd)
#
#del ccd1
#del ccd2
#del ratio
#
mask_as_ccd = CCDData(data=maskr.astype('uint8'), unit=u.dimensionless_unscaled)
mask_as_ccd.header['imagetyp'] = 'flat mask'
#
mask_as_ccd.write("../mask_new.fits")
