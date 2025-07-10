""" This script will clean the FFI images, it has been significantly reduced to only
show specific steps. The majority of the functional library can be found in the libraries
directory in utils, fits, photometry, image, and preprocessing libraries."""
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from config import Configuration
import os
from astropy.io import fits
import time
import numpy as np


class Clean:

    @staticmethod
    def clean_images():
        """ This is the main function script to clean multiple images, alternatively clean_img can be used to clean
        a single image.

        return no value is returned, the values images from in_path are cleaned and deposited in out_path
        """
        st = time.time()  # clock started

        # get the file list for all dates the FIELD was observed
        Utils.log("Getting file list...", "info")
        files, date_dirs = Utils.get_all_files_per_field(Configuration.RAW_DIRECTORY,
                                                         Configuration.FIELD,
                                                         'raw',
                                                         Configuration.FILE_EXTENSION)

        # make the output directories for the clean, difference, and flux files
        output_dirs = []

        for dte in date_dirs:
            output_dirs.append(Configuration.DATA_DIRECTORY + "clean/" + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + "clean/" + dte + "/" + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + "diff/" + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + "diff/" + dte + "/" + Configuration.FIELD)
            output_dirs.append(Configuration.DATA_DIRECTORY + "flux/" + dte)
            output_dirs.append(Configuration.DATA_DIRECTORY + "flux/" + dte + "/" + Configuration.FIELD)

        Utils.create_directories(output_dirs)
        # break if there are no files
        if len(files) == 0:
            Utils.log("No .fits files found for " + Configuration.FIELD + "!" +  ". Breaking...",
                      "debug")
            return()
        
        Utils.log("Starting to clean " + str(len(files)) + " images.", "info")
        for idx, file in enumerate(files):

            # make a new name for the file based on which actions are taken
            file_name = Preprocessing.mk_nme(file,
                                             'N',
                                             Configuration.SUBTRACT_BIAS,
                                             Configuration.SUBTRACT_DARK,
                                             Configuration.DIVIDE_FLAT,
                                             Configuration.CLIP_IMAGE,
                                             Configuration.SUBTRACT_SKY,
                                             Configuration.PLATE_SOLVE)

            # only create the files that don't exist
            if os.path.isfile(file_name) == 1:
                Utils.log("Image " + file_name +
                          " already exists. Skipping for now...", "info")

            # if the image does not exist then clean
            if os.path.isfile(file_name) == 0:

                # clean the image
                clean_img, header, bd_flag = Clean.clean_img(file)

                # write out the file
                if bd_flag == 0:
                    fits.writeto(file_name,
                                 clean_img, header, overwrite=True)

                    # print an update to the cleaning process
                    Utils.log("Cleaned image written as " + file_name + ".", "info")
                else:
                    Utils.log(file_name + " is a bad image. Not written.", "info")

            Utils.log(str(len(files) - idx - 1) + " images remain to be cleaned.",  "info")

        fn = time.time()  # clock stopped
        Utils.log("Imaging cleaning complete in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

    @staticmethod
    def clean_img(file):

        """ This function is the primary script to clean the image, various other functions found in this class
        can be found in the various libraries imported.

        :parameter  file - The file name of the image you would like to clean

        """

        Utils.log("Now cleaning " + file + ".", "info")

        # read in the image
        img, header = fits.getdata(file, header=True)

        # remove bias and dark as necessary
        if (Configuration.SUBTRACT_BIAS == 'Y') & (Configuration.SUBTRACT_DARK == 'Y'):
            st = time.time()
            bias, dark = Preprocessing.mk_combined_bias_and_dark(image_overwrite='N')
            img, header = Preprocessing.subtract_scaled_bias_dark(img, header)
            fn = time.time()
            Utils.log("Image bias and dark corrected in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping bias and darks subtraction...", "info")

        # flat divide if necessary
        if Configuration.DIVIDE_FLAT == 'Y':
            st = time.time()
            flat = Preprocessing.mk_flat(5, 300)
            img, header = Preprocessing.flat_divide(img, header)
            fn = time.time()
            Utils.log("Image flattened in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping image flattening....", "info")

        if Configuration.CLIP_IMAGE == 'Y':
            st = time.time()
            img, header = Preprocessing.clip_image(img, header)
            fn = time.time()
            Utils.log("Image over scan removed in " + str(np.around(fn - st, decimals=2)) + "s.", "info")
        else:
            Utils.log("Skipping over scan removal...", "info")

        # sky subtract if necessary
        if Configuration.SUBTRACT_SKY == 'Y':
            st = time.time()
            Utils.log("A background box of " + str(Configuration.PIX) + " x " + str(Configuration.PIX) +
                      " will be used for background subtraction.", "info")

            img, header = Preprocessing.sky_subtract(img, header, Configuration.WRITE_SKY)
            fn = time.time()
            Utils.log("Sky subtracted in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

        else:
            Utils.log("Skipping sky subtraction...", "info")

        if Configuration.PLATE_SOLVE == 'Y':
            st = time.time()
            Utils.log("Now plate solving and correcting the header.", "info")
            img, header = Preprocessing.correct_header(img, header)
            fn = time.time()
            Utils.log("The image has been plate sovled in " + str(np.around((fn - st), decimals=2)) + "s.", "info")
        Utils.log("Cleaning finished.", "info")

        bd_flag = 0

        return img, header, bd_flag
