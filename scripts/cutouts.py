from config import Configuration
from libraries.utils import Utils

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astroquery.vizier import Vizier
from ccdproc import ImageFileCollection
from matplotlib import pyplot as plt
from multiprocessing import Pool
from itertools import repeat
from tqdm import tqdm

import astropy.coordinates as coord
import astropy.units as u

import os

# The general logic for the main script
#   1. generate galaxy list
#   2. make master cutouts
#   3. make science cutouts
#   4. make differenced cutouts
#   5. save fits file with:
#       - WCS
#       - master cutout
#       - science cutout
#       - difference cutout

CUTOUT_WIDTH = 276 # in pixels
PARALLEL = False
PROCESSORS = Configuration.PROCESSORS

class Field:
    def __init__(self):
        self.name = Configuration.FIELD
        self.RA = Configuration.RA
        self.DEC = Configuration.DEC
        self.width = Configuration.FOV * 1.05 * u.deg
        self.coords = coord.SkyCoord(ra=self.RA, dec=self.DEC,
                                     unit=(u.deg, u.deg))

    def catalog_query(self):
        # TODO: check if local version exists otherwise
        # run query online
        vizier = Vizier(row_limit=-1) # minus 1 returns all rows
        results = vizier.query_region(self.coords,
                                      width=self.width,
                                      catalog='VII/291/gladep')
        return results[0] # only return the table

def single_cutout_from_filepath(filepath, position, output_dir, galaxy_id):
    """ Writes new cutout or stamp. Follow similar steps to
    https://docs.astropy.org/en/stable/nddata/utils.html#cutout-images

    :parameter filepath - path to file
    :parameter position - center position of the cutout
    :output_dir - 
    :galaxy_id - (str) galaxy id  

    :return Nothing is returned, the image is cutout

     """
    with fits.open(filepath) as hdul:
        hdu = hdul[0]
        wcs = WCS(hdu.header)

        # make the full path for cutout file
        cutout_path = output_dir
        cutout_path = os.path.dirname(cutout_path)
        datetimeinfo = hdu.header['DATE']
        datetimeinfo = datetimeinfo.replace("-", "")
        datetimeinfo = datetimeinfo.replace(":", "")
        cutout_name = "gp_"+galaxy_id+"_"+datetimeinfo+Configuration.FILE_EXTENSION
        cutout_filename = cutout_path + "/" + cutout_name
        Utils.log(f"Cutout filename: {cutout_filename}", "info")
        cutout_hdu = fits.PrimaryHDU(data=hdu.data, header=hdu.header)
        if not os.path.isfile(cutout_filename):
            try:
                cutout = Cutout2D(hdu.data, position, CUTOUT_WIDTH, wcs=wcs, mode='partial')
                # save cutout
                cutout_hdu.data = cutout.data
                cutout_hdu.header.update(cutout.wcs.to_header())
                cutout_hdu.writeto(cutout_filename)
            except NoOverlapError:
                Utils.log(f"No Overlap!! Skipping.", "info")
        else:
            Utils.log(f"Cutout {cutout_filename} already exists. Skipping...", "info")



def generate_master_cutouts():
    # get list of galaxies based on toros field from GLADE+
    toros_field = Field()
    galaxies = toros_field.catalog_query()

    # consider galaxies on edges of frame

    # load master frame
    file = Configuration.MASTER_DIRECTORY + toros_field.name + "_master" + Configuration.FILE_EXTENSION
    hdu = fits.open(file)[0]
    wcs = WCS(hdu.header)

    if not wcs.is_celestial:
        wcs = None

    for galaxy in galaxies:
        position = coord.SkyCoord(ra=galaxy['RAJ2000']*u.deg, dec=galaxy['DEJ2000']*u.deg)
        galaxy_id = galaxy["GLADE+"].astype("str").zfill(9)
        single_cutout_from_filepath(file, position, Configuration.CUTOUT_DIRECTORY+toros_field.name+ "/", galaxy_id)

def generate_cutouts_from_filepath(filepath, filetype=""):
    Utils.log(f"Generating {filetype} cutouts for file: {filepath}", "info")
    # get list of galaxies based on toros field from GLADE+
    toros_field = Field()
    galaxies = toros_field.catalog_query()

    # consider galaxies on edges of frame

    # make output directroy
    output_dir = Configuration.CUTOUT_DIRECTORY + toros_field.name + "/"
    if filetype.lower() == "diff":
        output_dir += "diff/"
    elif filetype.lower() == "master":
        output_dir += "master/"
    elif filetype.lower() == "clean":
        output_dir += "clean/"
    else:
        Utils.log(f"Filetype {filetype.lower()} is not implemented. Writing files to {output_dir}.", "info")

    # load master frame
    hdu = fits.open(filepath)[0]
    wcs = WCS(hdu.header)

    if not wcs.is_celestial:
        wcs = None

    if PARALLEL:
        Utils.log("Using multiprocessing!","info")
        # Generate iterable of positions for each galaxy
        positions = []
        galaxy_ids = []
        for galaxy in galaxies:
            position = coord.SkyCoord(ra=galaxy['RAJ2000']*u.deg,
                                      dec=galaxy['DEJ2000']*u.deg)
            positions.append(position)
            galaxy_ids.append(galaxy["GLADE+"].astype("str").zfill(9))
        # pack args for starmap
        args = zip(repeat(filepath), positions, repeat(output_dir), galaxy_ids)

        with Pool(processes=PROCESSORS) as pool:
            pool.starmap(single_cutout_from_filepath, args)
    else:
        for galaxy in galaxies:
            position = coord.SkyCoord(ra=galaxy['RAJ2000']*u.deg,
                                      dec=galaxy['DEJ2000']*u.deg)
            galaxy_id = galaxy["GLADE+"].astype("str").zfill(9)
            single_cutout_from_filepath(filepath,  position, output_dir, galaxy_id)

    Utils.log(f"Done!! Cutout process complete for {filepath}", "info")

def fits_to_jpg(filename):
    # Load data
    data = fits.getdata(filename)
    zscale = ZScaleInterval()
    plt.imshow(zscale(data), cmap='gray')
    filename = os.path.splitext(filename)[0] + ".png"
    plt.axis('off')
    plt.savefig(filename, bbox_inches='tight', pad_inches=0)

def generate_cutouts_clean_per_field():
    fieldname = Field().name
    Utils.log(f"Generating cutouts for all CLEAN files for {fieldname}", "info")
    files, dates = Utils.get_all_files_per_field(
                    Configuration.CLEAN_DIRECTORY,
                    fieldname,
                    'cutout',
                    '.fits')

    numfiles = len(files)

    for idx, file in enumerate(tqdm(files)):
        generate_cutouts_from_filepath(file, filetype="clean")
        Utils.log(f"{idx} of {numfiles} files complete","info")
    Utils.log("All done making CLEAN cutouts!", "info")

def generate_all_cutouts_per_field(filetype="clean"):
    fieldname = Field().name
    Utils.log(f"Generating cutouts for all {filetype.upper()} files for {fieldname}", "info")
    if filetype.lower() == "clean":
        directory = Configuration.CLEAN_DIRECTORY
    elif filetype.lower() == "diff":
        directory = Configuration.DIFFERENCED_DIRECTORY
    else:
        Utils.log(f"Filetype {filetype} is not known.", "info")

    files, dates = Utils.get_all_files_per_field(
                    directory,
                    fieldname,
                    'cutout',
                    '.fits')

    numfiles = len(files)

    for idx, file in enumerate(tqdm(files)):
        generate_cutouts_from_filepath(file, filetype="clean")
        Utils.log(f"{idx} of {numfiles} files complete","info")
    Utils.log("All done making CLEAN cutouts!", "info")


if __name__ == "__main__":
    pass
