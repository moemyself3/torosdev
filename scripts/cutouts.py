from config import Configuration
from libraries.utils import Utils

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from ccdproc import ImageFileCollection

import astropy.coordinates as coord
import astropy.units as u

import os

CUTOUT_WIDTH = 276 # in pixels
PARALLEL = False

class Field:
    def __init__(self):
        self.name = Configuration.FIELD
        self.RA = Configuration.RA
        self.DEC = Configuration.DEC
        self.width = Configuration.FOV * 1.05 * u.deg
        self.coords = coord.SkyCoord(ra=self.RA, dec=self.DEC,
                                     unit=(u.deg, u.deg))

    def catalog_query(self):
        vizier = Vizier()
        results = vizier.query_region(self.coords,
                                      width=self.width,
                                      catalog='VII/291/gladep')
        return results[0] # only return the table

def single_cutout(hdu, wcs, position, differenced_path, galaxy):
    """ Writes new cutout or stamp. Follow similar steps to
    https://docs.astropy.org/en/stable/nddata/utils.html#cutout-images

    :parameter hdu - The hdu to make the cutout from. Contains data & header
    :parameter wcs - WCS object of the original frame.
    :parameter position - center position of the cutout

    :return Nothing is returned, the image is cutout

     """
    # make the full path for cutout file
    cutout_path = differenced_path.replace("/diff/", "/cutout/")
    cutout_path = os.path.dirname(cutout_path)

    datetimeinfo = hdu.header['DATE']
    cutout_name = "gp_"+galaxy["GLADE+"].zfill(9)+"_"+datetimeinfo+Configuration.FILE_EXTENSION
    cutout_filename = cutout_path + "/" + cutout_name

    if not os.isfile(cutout_filename):
        cutout = Cutout2D(hdu.data, position, CUTOUT_WIDTH, wcs=wcs)
        # save cutout
        cutout.write(cutout_path)
    else:
        Utils.log(f"Cutout {cutout_path} already exists. Skipping...", "info")

def generate_cutouts():
    # get list of galaxies based on toros field from GLADE+
    toros_field = Field()
    galaxies = toros_field.catalog_query()

    # consider galaxies on edges of frame

    # load master, source, and differenced frames
    hdu = fits.open(file)[0]
    wcs = WCS(hdu.header)

    if not wcs.is_celestial:
        wcs = None

    if PARALLEL:
        # pack args for starmap
        args = None
        with Pool(processes=4) as pool:
            pool.starmap(cutout, args)
    else:
        for galaxy in galaxies:
            cutout(hdu, wcs, postion, cutout_name)


if __name__ == "__main__":
    pass
