from config import Configuration
from libraries.utils import Utils

from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS as pywcs
from astropy.io import fits
from ccdproc import ImageFileCollection

clean_files, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY, Configuration.FIELD, "clean", Configuration.FILE_EXTENSION)
collection = ImageFileCollection(filenames=clean_files)

bad_list = []
for header, file_name in collection.headers(return_fname=True):
    wcs = pywcs(header)
    x_scale, y_scale = proj_plane_pixel_scales(wcs) * 60 * 60
    if x_scale > 0.5 or y_scale > 0.5:
        bad_list.append(file_name)
        print(f"file_name: {file_name}")
        print(f"date: {header['date']}")
        print(f"platescale: {x_scale}, {y_scale}")

