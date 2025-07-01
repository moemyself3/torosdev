import pandas as pd
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
from libraries.utils import Utils

image, header = fits.getdata('/Users/yuw816/Data/toros/commissioning/clean/2024-11-11/FIELD_0e.001/FIELD_0e.001_300s_12000x10600_1239_bkfcsp.fits', header=True)
erosita = pd.read_csv("~/Development/toros/analysis/erosita_sources.csv", sep=',')

erosita['ra_deg'] = erosita.apply(lambda x: (x.rah + x.ram / 60 + x.ras / 3600) * 360 / 24, axis=1)
erosita['de_deg'] = erosita.apply(lambda x: x.ded - x.dem / 60  - x.des / 3600, axis=1)

# convert the ra & dec positions to x & y positions based on the image header
w = WCS(header)
ra = erosita.ra_deg.to_numpy()
dec = erosita.de_deg.to_numpy()

# convert to x, y
x, y = w.all_world2pix(ra, dec, 0)
positions = (x, y)

# add the x/y to the star data frame
erosita['x'] = x
erosita['y'] = y

f = open("/Users/yuw816/Development/toros/analysis/erosita.reg", "w")
for xx, yy in zip(x, y):
    f.write("circle(" + str(np.around(xx, decimals=2)) + ',' + str(np.around(yy, decimals=2))+", 30)\n")
f.close()
