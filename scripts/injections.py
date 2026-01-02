""" This script will create training set for ML by injecting sources using the 
ePSF of the image supplied.
"""

from libraries.utils import Utils
from libraries.injector import Injector
from scripts.difference import BigDiff
from config import Configuration

from astropy.io import fits
from astropy.table import QTable

import numpy as np
import pandas as pd

import os
import subprocess
import time

def create_injection_images() -> None:
    """ This function creates the ePSF model for images and injects sources
    at random positions with at random magnitudes based on the range of 
    magnitudes of the sources found in the image.
    """

    # get the file list for all dates the FIELD was observed
    Utils.log("Creating injection images... Getting file list...", "info")
    files, date_dirs = Utils.get_all_files_per_field(
            Configuration.CLEAN_DIRECTORY,
            Configuration.FIELD,
            'inject',
            Configuration.FILE_EXTENSION)

    # make the output directories
    output_dirs = []
    for date in date_dirs:
        output_dirs.append(Configuration.EPSF_DIRECTORY)
        output_dirs.append(Configuration.EPSF_MODELS_DIRECTORY)
        output_dirs.append(Configuration.EPSF_CANDIDATES_DIRECTORY)
        output_dirs.append(Configuration.EPSF_MODELS_DIRECTORY      + date)
        output_dirs.append(Configuration.EPSF_CANDIDATES_DIRECTORY  + date)
        output_dirs.append(Configuration.INJECTED_CLEAN_DIRECTORY   + date)
        output_dirs.append(Configuration.INJECTED_DIFF_DIRECTORY    + date)
        output_dirs.append(Configuration.INJECTED_FLUX_DIRECTORY    + date)
        output_dirs.append(Configuration.INJECTED_CATALOG_DIRECTORY + date)
        output_dirs.append(
                Configuration.EPSF_MODELS_DIRECTORY     + date \
                        + "/" + Configuration.FIELD)
        output_dirs.append(
                Configuration.EPSF_CANDIDATES_DIRECTORY + date \
                        + "/" + Configuration.FIELD)
        output_dirs.append(
                Configuration.INJECTED_CLEAN_DIRECTORY  + date \
                        + "/" + Configuration.FIELD)
        output_dirs.append(
                Configuration.INJECTED_DIFF_DIRECTORY   + date \
                        + "/" + Configuration.FIELD)
        output_dirs.append(
                Configuration.INJECTED_FLUX_DIRECTORY   + date \
                        + "/" + Configuration.FIELD)
        output_dirs.append(
                Configuration.INJECTED_CATALOG_DIRECTORY+ date \
                        + "/" + Configuration.FIELD)

    Utils.create_directories(output_dirs)
    # break if there are no files
    if len(files) == 0:
        Utils.log(
                "No .fits files found for " + Configuration.FIELD + "!" \
                        +  ". Breaking...",
                  "debug")
        return()

    # use files list and update paths to training/injected/clean as output
    # run loop for each file to inject sources and save new injected clean frame
    for file in files:
        Utils.log(f"injection process started on {file}...", "info")

        injected_filepath = file.replace(
                                Configuration.CLEAN_DIRECTORY,
                                Configuration.INJECTED_CLEAN_DIRECTORY)
        injected_dir, injected_filename = os.path.split(injected_filepath)
        injected_filename = os.path.splitext(injected_filename)[0] + 'i' \
                            + Configuration.FILE_EXTENSION
        injected_filepath = injected_dir + os.path.sep \
                            + injected_filename

        if os.path.isfile(injected_filepath):
            Utils.log(
                    f"Skipping... Injected file {injected_filepath} exists",
                    "info")
            continue

        Utils.log(f"working on injections for {file}", "info")

        # get file name differnced image catalog
        diff_catalog = file.replace(
                            Configuration.CLEAN_DIRECTORY,
                            Configuration.DIFFERENCED_CATALOG_DIRECTORY)
        diff_catalog = os.path.splitext(diff_catalog)[0] + "ad.cat"
        Utils.log(f"Reading {diff_catalog}","info")

        # check that catalog exists
        if not os.path.isfile(diff_catalog):
            Utils.log("Skipping. No catalog found.", "info")
            continue # skip to next file

        # get number of sources from differenced image
        diff_sources = make_source_extractor_catalog_dataframe(diff_catalog)
        num_sources = len(diff_sources)

        # make sure there are sources in catalog
        if num_sources == 0:
            Utils.log(
                    "Skipping injections. No sources in the catalog...", 
                    "info")
            continue # skip to the next file

        # start injector obect and make epsf filepath
        injector = Injector(file)
        epsf_path = make_epsf_filepath(file)
        candidate_starlist_path = make_candidate_starlist_filepath(file)
        # check if candidate starlist exists
        if not os.path.isfile(candidate_starlist_path):
            Utils.log("Making candidate starlist.", "info")
            lonely_peaks_table = injector.make_candidate_starlist() # takes ~4 mins
            save_candidate_starlist(lonely_peaks_table, file)
        else:
            lonely_peaks_table = load_candidate_starlist(candidate_starlist_path)
            # Update Injector with lonely_peaks_table
            injector.lonely_peaks_table = lonely_peaks_table
            injector.image_median = lonely_peaks_table.meta['image_median']

        # check if EPSF file exists
        if not os.path.isfile(epsf_path):
            Utils.log("EPSF model not found. Making new model.", "info")
            injector.make_epsf()
            Utils.log("epsf computed", "info")
            save_epsf(injector, file)
        else:
            Utils.log("EPSF model found. Using existing model.", "info")
            # load epsf model called epsf_resized from Injector
            epsf_resized = fits.getdata(epsf_path)
            injector.load_epsf(epsf_resized)

        # get magnitude range from loneley peaks and add +/- 0.5
        upper_magnitude = lonely_peaks_table['mag'].min() - 0.5
        lower_magnitude = lonely_peaks_table['mag'].max() + 0.5
        upper_mag_limit = Configuration.UPPER_MAG_LIMIT - Configuration.MAG_OFFSET
        # REMEMBER the lower the magnitude the BRIGHTER the object
        # check if the upper_magnitude of lonely_peaks_table is BRIGHTER
        #   than the upper_mag_limit
        if upper_magnitude < upper_mag_limit:
            upper_magnitude = upper_mag_limit

        # make path and name for injection catalog
        injection_catalog = file.replace(
                                Configuration.CLEAN_DIRECTORY,
                                Configuration.INJECTED_CATALOG_DIRECTORY)
        injection_catalog = injection_catalog.replace(".fits", "_injection.csv")

        # Check for injection catalog
        if not os.path.isfile(injection_catalog):
            # Generate N random sources to match with sources found in 
            Utils.log(
                    f"Making {num_sources} sources using upper mag {upper_magnitude} and lower mag {lower_magnitude}",
                    "info")
            positions = generate_random_source_positions(num_sources)
            position_table = pd.DataFrame(positions, columns=['X','Y'])
            position_table['MAG'] = generate_random_source_magnitudes(
                                        upper_magnitude,
                                        lower_magnitude,
                                        num_sources)

            # save injection source table
            position_table.to_csv(injection_catalog, index=False, sep="\t")
        else:
            Utils.log(f"Using existing injection catalog.", "info")
            position_table = pd.read_csv(injection_catalog, sep="\t")

        # inject sources using position table X, Y, MAG
        for index, source in position_table.iterrows():
            position = (source['X'], source['Y'])
            injector.inject(position, source['MAG'], injector.image.data)

        # Add injection information to header
        injector.image.header.append(
                ('INJECTED', num_sources, 'number of sources injected')
                )
        Utils.log(f"Writing injected file to {injected_filepath}", "info")
        # Save injected fits file to injected CLEAN directory
        injector.image.writeto(injected_filepath)
        injector.hdul.close()

def make_candidate_starlist_filepath(file:str) -> str:
    filename = os.path.basename(file)
    starlist_path = file.replace(Configuration.CLEAN_DIRECTORY,
                             Configuration.EPSF_CANDIDATES_DIRECTORY)
    starlist_path = starlist_path.replace(filename, "STARS_"+filename)
    ext = ".ecsv"
    starlist_path = os.path.splitext(starlist_path)[0] + ext
    return starlist_path

def save_candidate_starlist(candidate_starlist: QTable, file) -> None:
    starlist_path = make_candidate_starlist_filepath(file) 
    candidate_starlist.write(starlist_path, format='ascii.ecsv')
    Utils.log(f"EPSF star candidates table saved {starlist_path}", "info")

def load_candidate_starlist(candidate_starlist_path: str) -> QTable:
    Utils.log("Loading EPSF candidate starlist.", "info")
    return QTable.read(candidate_starlist_path, format="ascii.ecsv")

def make_epsf_filepath(file:str) -> str:
    filename = os.path.basename(file)
    epsf_path = file.replace(Configuration.CLEAN_DIRECTORY,
                             Configuration.EPSF_MODELS_DIRECTORY)
    epsf_path = epsf_path.replace(filename, "EPSF_"+filename)
    return epsf_path

def save_epsf(injector, file) -> None:
    filename = os.path.basename(file)
    epsf_path = make_epsf_filepath(file)
    epsf_filename = os.path.basename(epsf_path)
    epsf_model = injector.epsf_resized
    hdu = fits.PrimaryHDU(data=epsf_model, header=None)
    nstars = len(injector.stars.all_good_stars)
    hdu.header.append(
            ('PSFSTARS', nstars, 'Number of stars used for psf fitting' ))
    hdu.header.append(
            ('SOURCEIM', filename, 'source image'))
    fits.writeto(epsf_path, hdu.data, hdu.header, overwrite=True)
    Utils.log(f"EPSF model saved for {file}", "info")

def make_difference_image_catalog() -> None:
    # get the file list for all dates the FIELD was observed
    Utils.log("Getting diff file list...", "info")
    files, date_dirs = Utils.get_all_files_per_field(
            Configuration.DIFFERENCED_DIRECTORY,
            Configuration.FIELD,
            'inject',
            Configuration.FILE_EXTENSION)

    # make the output directories for the differenced catalogs
    output_dirs = []
    for date in date_dirs:
        output_dirs.append(Configuration.DIFFERENCED_CATALOG_DIRECTORY + date )
        output_dirs.append(
                Configuration.DIFFERENCED_CATALOG_DIRECTORY \
                        + date + "/"+ Configuration.FIELD)

    Utils.create_directories(output_dirs)
    Utils.log(
            f"Making differenced image catalogs for {Configuration.FIELD}",
            "info")
    for file in files:
        catalog_file = file.replace("/diff/", "/diff_catalogs/")
        # get file name without extension
        catalog_file_name = os.path.splitext(catalog_file)[0]
        catalog_file_name = catalog_file_name + ".cat"
        catalog_file = catalog_file_name

        if not os.path.exists(catalog_file):
             # run source extractor to make source catalogs
             source_extractor_config = os.getcwd()+ '/' + "rf.sex"
             command = Configuration.SOURCE_EXTRACTOR_PATH + " " + file \
                         + " -c " + source_extractor_config \
                         + " -CATALOG_NAME " + catalog_file
             subprocess.run(command, shell=True)
        else:
            Utils.log(f"Skipping... Catalog {catalog_file_name} exists!","info")

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

def generate_random_source_positions(
        num_sources: int
) -> list[list[float, float]]:
    rng = np.random.default_rng() # Random Number Generator
    # uniform distribution of num_sources positions
    positions = rng.uniform(31, Configuration.AXS_X-31, (num_sources, 2))
    return positions

def generate_random_source_magnitudes(
        upper_magnitude: float,
        lower_magnitude: float,
        num_sources: int
) -> list[float]:
    rng = np.random.default_rng() # Random Number Generator
    # uniform distribution of magnitudes
    magnitudes = rng.uniform(upper_magnitude, lower_magnitude, num_sources)
    return magnitudes

def update_config_path_to_injected_dir() -> None:
    Utils.log("Updating DATA_DIRECTORY and CLEAN_DIRECTORY.", "info")
    Configuration.DATA_DIRECTORY = Configuration.INJECTED_DATA_DIRECTORY
    Configuration.CLEAN_DIRECTORY = Configuration.DATA_DIRECTORY + "clean/"
    Configuration.DIFFERENCED_DIRECTORY = Configuration.DATA_DIRECTORY + "diff/"

def difference_injected_images() -> None:
    # Similar process used in difference.py difference_images
    # This uses injected clean images. 
    # difference_images uses the default clean dir. need to use injected/clean
    Utils.log("Differencing injeced images...", "info")
    BigDiff.difference_images()

def main() -> None:
    make_difference_image_catalog()
    create_injection_images()
    update_config_path_to_injected_dir()
    difference_injected_images()


if __name__ == "__main__":
    main()
