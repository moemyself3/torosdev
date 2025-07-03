from libraries.utils import Utils
from scripts.lightcurves import Lightcurves
from scripts.clean import Clean
from scripts.master import Master
from scripts.difference import BigDiff
from config import Configuration

# do necessary prep work such as making output directories
Utils.create_directories(Configuration.DIRECTORIES)

# do the necessary preprocessing of the images
if Configuration.CLEAN_SKIP == 'N':
    Clean.clean_images(image_clip='Y', bias_subtract='Y', dark_subtract='Y',
                       flat_divide='Y', sky_subtract='Y', plate_solve='Y')
else:
    Utils.log("Skipping image cleaning.", "info")

if Configuration.MASTER_SKIP == 'N':
    master, star_list = Master.pull_master()
else:
    Utils.log("Skipping master frame generation.", "info")

if Configuration.DIFFERENCE_SKIP == 'N':
    BigDiff.difference_images()
else:
    Utils.log("Skipping image differencing.", "info")

if Configuration.PHOTOMETRY_SKIP == 'N':
    Lightcurves.generate_flux_files()
else:
    Utils.log("Skipping photometry.", "info")

if Configuration.LIGHTCURVE_SKIP == 'N':
    Lightcurves.mk_raw_lightcurves(star_list)
else:
    Utils.log("Skipping making raw light curves.", "info")
Utils.log("All done! See ya later, alligator.", "info")
