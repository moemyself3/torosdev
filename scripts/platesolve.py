from config import Configuration
from libraries.preprocessing import Preprocessing
from libraries.utils import Utils
from astropy.io import fits

import os
import subprocess
import shutil

def get_files():
    # get the file list for all dates the FIELD was observed
    Utils.log("Getting file list...", "info")
    files, date_dirs = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY,
                                                     Configuration.FIELD,
                                                     'raw',
                                                     Configuration.FILE_EXTENSION)

    # break if there are no files
    if len(files) == 0:
        Utils.log("No .fits files found for " + Configuration.FIELD + "!" +  ". Breaking...",
                  "debug")
    return files


def process_files(files):
    Utils.log("Starting to process " + str(len(files)) + " images.", "info")
    for idx, file in enumerate(files):
         # read in the image
        img, header = fits.getdata(file, header=True)

        # make a new name for the file based on which actions are taken
        file_name = Preprocessing.mk_nme(file,
                                         'N',
                                         Configuration.SUBTRACT_BIAS,
                                         Configuration.SUBTRACT_DARK,
                                         Configuration.DIVIDE_FLAT,
                                         Configuration.CLIP_IMAGE,
                                         Configuration.SUBTRACT_SKY,
                                         Configuration.PLATE_SOLVE)
        file_temp = file_name.split("_")
        file_temp.pop(-2)
        file_name = "_".join(file_temp)
        Utils.log("Now plate solving and correcting the header.", "info")
        img, header = correct_header(img, header, file_name)
        Utils.log(str(len(files) - idx - 1) + " images remain to be cleaned.",  "info")



def correct_header(img, header, clean_file_name=None):
    if Configuration.PLATE_SOLVE_METHOD.lower() == 'astrometry_net' and clean_file_name is not None:
        Utils.log("Plate solve method: astrometry_net", "info")

        working_directory_file = Configuration.WORKING_DIRECTORY + "analysis/" + os.path.basename(clean_file_name) +".temp"
 
        fits.writeto(working_directory_file, img, header, overwrite=True)
        output_file = os.path.splitext(working_directory_file)[0]

        # set output directory for astrometry
        output_dir = os.path.dirname(working_directory_file)

        astrometry_command = "solve-field --use-source-extractor --source-extractor-path "+ Configuration.SOURCE_EXTRACTOR_PATH + " --scale-units arcsecperpix --scale-low 0.48 --scale-high 0.5 --no-plots --temp-axy --index-xyls none --match none --rdls none --solved none --corr none --dir " + output_dir + " --new-fits " + output_file + " " + working_directory_file

        subprocess.run(astrometry_command, shell=True)

        # Clean up astrometry output, remove .temp and .wcs files
        astrometry_extensions = [ ".temp", ".wcs"]
        
        # First check that astrometry generated plate solved file (.fits)
        if os.path.exists(output_file):
            shutil.move(output_file, clean_file_name)
            Utils.log(f"Cleaned image written as " + clean_file_name + ".", "info")
            # if file exists then remove extra files
            for extension in astrometry_extensions:
                file_to_remove = output_file + extension
                if os.path.exists(file_to_remove):
                    os.remove(file_to_remove)
                    Utils.log(f"File {file_to_remove} deleted.", "info")
        else:
            # if file does not exist then move the temp file out of clean for review
            temp_file_name = clean_file_name.replace("/clean/","/review/")
            shutil.move(working_directory_file, temp_file_name)
            Utils.log(f"Plate Solve Failed. Moving file to review directory for human inspection.", "info") 
    else:
        Utils.log("not using astrometry_net...", "info")

    return img, header


if __name__ == "__main__":
    files = get_files()
    process_files(files)
