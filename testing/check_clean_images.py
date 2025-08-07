from config import Configuration
from libraries.utils import Utils

from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS as pywcs
from astropy.io import fits
from ccdproc import ImageFileCollection

import os

clean_files, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY, Configuration.FIELD, "clean", Configuration.FILE_EXTENSION)
collection = ImageFileCollection(filenames=clean_files)

def check_platescale():
    bad_list = []
    for header, file_name in collection.headers(return_fname=True):
        wcs = pywcs(header)
        x_scale, y_scale = proj_plane_pixel_scales(wcs) * 60 * 60
        if x_scale > 0.5 or y_scale > 0.5:
            bad_list.append(file_name)
            #print(f"file_name: {file_name}")
            #print(f"date: {header['date']}")
            #print(f"platescale: {x_scale}, {y_scale}")
    print(f"total images: {len(collection.files)}")
    print(f"total bad wcs: {len(bad_list)}")
    return bad_list

def check_wcs_exists():
    return collection.files_filtered(wcsaxes=None, include_path=True)

def check_bad_wcs():
    return collection.files_filtered(bad_wcs='Y')

def check_badimage():
    return collection.files_filtered(badimage='Y')

def rm_badscale():
    good_wcs_path = Configuration.CLEAN_DIRECTORY + "2024-11-02/" + Configuration.FIELD + "/"
    good_wcs_files = os.listdir(good_wcs_path)
    good_wcs_files.sort()
    good_wcs_path = good_wcs_path + good_wcs_files[0]

    goodwcs = fits.getheader(good_wcs_path)
    wcs = pywcs(goodwcs)

    badscalepath = []
    bad_scale = check_platescale()
    for bad_file in bad_scale:
        for clean_file in clean_files:
            if bad_file in clean_file:
                badscalepath.append(clean_file)

    for badwcs in badscalepath:
        with fits.open(badwcs, mode="update") as hdul:
            hdr = hdul[0].header
            for item in wcs.to_header():
                try:
                    del hdr[item]
                except KeyError:
                    pass
            hdul.flush()
