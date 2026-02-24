from config import Configuration
from libraries.utils import Utils

from astroquery.vizier import Vizier
from astropy.wcs import WCS
from astropy.io import fits

import astropy.coordinates as coord
import astropy.units as u
import pandas as pd
import numpy as np
import os

CATALOGS = [
        'B/gcvs/gcvs_cat',  # GENERAL CATALOG OF VARIABLE STARS
        'B/vsx/vsx',        # INTERNATIONAL VARIABLE STAR INDEX
        ]

def generate_directories():
    # get the file list for all dates the FIELD was observed
    Utils.log("Generating directories for Real Bogus training set", "info")
    files, date_dirs = Utils.get_all_files_per_field(
            Configuration.INJECTED_CLEAN_DIRECTORY,
            Configuration.FIELD,
            'realbogus',
            Configuration.FILE_EXTENSION)

    # make the output directories
    output_dirs = []
    for date in date_dirs:
        output_dirs.append(Configuration.REALBOGUS_CATALOG_DIRECTORY)
        output_dirs.append(Configuration.REALBOGUS_CATALOG_DIRECTORY + date)
        output_dirs.append(
                Configuration.REALBOGUS_CATALOG_DIRECTORY   + date \
                        + "/" + Configuration.FIELD)

    Utils.create_directories(output_dirs)

    return files, date_dirs

def make_source_extractor_catalog_dataframe(file: str) -> pd.DataFrame:
   # get header information from the first set of rows
   # look '#' to find the header keys and
   # stop when the line does not start with '#'
    header = []
    with open(file, "r") as lines:
        for line in lines:
            if line.startswith("#"):
                header.append(line)
            else:
                break

    # reshape header keys to columns
    header_keys = [line.split()[2] for line in header]
    df = pd.read_csv(file, sep='\s+', skiprows=len(header_keys),
                     names=header_keys, header=None)
    return df

def make_labeled_catalog(cleanfilename, gcvs_coords, vsx_coords):
# Catalog labeling
# Cross match catalog with:
#   Variable stars as REAL
#   XY position of injection as REAL

    # load sources from differenced frame as dataframe
    diff_sources = load_diff_catalog(cleanfilename)

    # add real bogus flag real = 1,  bogus = 0
    #   use bogus as default
    diff_sources['REAL'] = 0

    # convert differenced sources RA DEC to skycoords
    diff_coords = format_source_coords(diff_sources)

    # crossmatch variables with differenced sources
    idx_matches = crossmatch_variable_star_catalogs(
            gcvs_coords,
            vsx_coords,
            diff_coords)

    # convert injection XY positions to RA DEC SkyCoords
    inj_coords = convert_injection_XY_to_world(cleanfilename)

    # crossmatch injections with differenced sources
    inj_matches = crossmatch_catalogs(inj_coords, diff_coords)

    # join IDs for Variable Stars matches and Injection matches
    idx_matches = np.concatenate((idx_matches, inj_matches), axis=0)

    # make unique list of ids
    idx_matches = list(set(idx_matches))

    # add real flag to data
    diff_sources.loc[idx_matches, 'REAL'] = 1

    # save file for training
    Utils.log(f"saving file: {realbogusfilename}", "info")
    diff_sources.to_csv(realbogusfilename, index=False)

def get_variable_star_table():
# Query FIELD for variable stars
#   results[0] is gcvs catalog
#   results[1] is vsx catalog
    ra = Configuration.RA * u.deg
    dec = Configuration.DEC * u.deg
    coords = coord.SkyCoord(ra=ra, dec=dec,
                            unit=(u.deg, u.deg))
    width = Configuration.FOV * 1.05 * u.deg
    vizier = Vizier(row_limit=-1) # minus 1 returns all rows
    variable_star_table = vizier.query_region(coords,
                                width=width,
                                catalog=CATALOGS)
    return variable_star_table

def format_gcvs_coords(variable_star_table):
    ra  = variable_star_table[0]["RAJ2000"]
    dec = variable_star_table[0]["DEJ2000"]
    return coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg))


def format_vsx_coords(variable_star_table):
    ra  = variable_star_table[1]["RAJ2000"]
    dec = variable_star_table[1]["DEJ2000"]
    return coord.SkyCoord(ra, dec, unit=(u.deg, u.deg))

def format_source_coords(dataframe):
    ra = dataframe['ALPHA_J2000']
    dec = dataframe['DELTA_J2000']
    return coord.SkyCoord(ra, dec, unit=(u.deg, u.deg))

def load_injection_catalog(filename):
# Load matching injection catalog for a given file.
# file name shares root up to bkfcsp 
#  eg. FIELD_300s_12000x10600_888_bkfcsp_injection.csv
    injection_catalog = make_source_catalog_name(filename, "injection")
    return pd.read_csv(injection_catalog, sep="\t")

def load_diff_catalog(filename):
    diff_catalog = make_source_catalog_name(
            filename, 
            source_catalog_type="difference")
    return  make_source_extractor_catalog_dataframe(diff_catalog)

def make_source_catalog_name(cleanfilename, source_catalog_type="injection"):
    cleanfile_ext = "i.fits"
    cleanfile_dir = "/clean/"
    source_catalogs = {
            "injection":
                {
                    "ext":  "_injection.csv",
                    "dir":  "/catalogs/",
                    },
            "difference":
                {
                    "ext":  "iad.cat",
                    "dir":  "/diff_catalogs/",
                    },
            "realbogus":
                {
                    "ext": "iad_realbogus.csv",
                    "dir": "/rb_catalogs/",
                    },
                }

    # change file extension
    filename = cleanfilename.replace(
            cleanfile_ext,
            source_catalogs[source_catalog_type]['ext']
            )
    # change directory
    filename = filename.replace(
            cleanfile_dir,
            source_catalogs[source_catalog_type]['dir'])

    return filename

def convert_injection_XY_to_world(filename):
    header = fits.getheader(filename)
    wcs = WCS(header=header)
    xy_catalog = load_injection_catalog(filename)
    x = xy_catalog['X'].to_numpy()
    y = xy_catalog['Y'].to_numpy()
    return wcs.pixel_to_world(x, y)

def crossmatch_variable_star_catalogs(gcvs_coords, vsx_coords, source_coords):
    gcvs_matches = crossmatch_catalogs(gcvs_coords, source_coords)
    vsx_matches  = crossmatch_catalogs(vsx_coords, source_coords)

    # make list of IDs between idx_gcvs and idx_vsx
    idx_matches = np.concatenate((gcvs_matches, vsx_matches), axis=0)

    # remove repeated idxs
    idx_matches = set(idx_matches) # makes a unique set of values
    idx_matches = list(idx_matches) # change set to list

    return idx_matches

def crossmatch_catalogs(catalog_coords, source_coords):
    idx, d2d, d3d = catalog_coords.match_to_catalog_sky(source_coords)

    # Apply an angular distance filter
    max_sep = 5.0 * u.arcsec

    return idx[d2d < max_sep]

def main():
    # setup directories for Real Bogus catalogs and get cleanfile list
    cleanfiles, clean_date_dirs = generate_directories()

    # get table of variable stars
    variable_star_table = get_variable_star_table()

    # format startable to get skycoords
    gcvs_coords = format_gcvs_coords(variable_star_table)
    vsx_coords = format_vsx_coords(variable_star_table)

    for cleanfilename in cleanfiles:
        Utils.log(f"working on {cleanfilename}", "info")
        realbogusfilename = make_source_catalog_name(cleanfilename, "realbogus")
        if os.path.isfile(realbogusfilename):
            Utils.log(f"rb_catalog already exists! Skipping...", "info")
            continue

        Utils.log(f'{cleanfilename=}',"info")
        make_labeled_catalog(cleanfilename, gcvs_coords, vsx_coords)
