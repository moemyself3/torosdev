""" This serves as the configuration file for the TOROS pipeline. """


class Configuration:

    # Computer for reduction
    MACHINE = 'm1'
    RAW_FILE_EXTENSION = '.fits'
    FILE_EXTENSION = '.fits'

    # update for different data products
    FIELD = 'FIELD_2b.022'
    RA = 52.759
    DEC = -35.611

    # is there a transient you want the light curve for?
    TRANSIENT_LC = 'N'
    TRANSIENT_NAME = 'AT2024xsq' # 'AT2024xhm'  # 'AT2024ykx' #
    TRANSIENT_RA = 10.1958012 # 46.9968525  # 307.422951333 #
    TRANSIENT_DEC = -21.9258309  # -38.93038375  # -30.269894 #

    # steps to skip
    CLEAN_SKIP = 'N'
    WRITE_SKY = 'Y'
    MASTER_SKIP = 'Y'
    DIFFERENCE_SKIP = 'Y'
    PHOTOMETRY_SKIP = 'Y'
    LIGHTCURVE_SKIP = 'Y'

    # how do you want to clean the image?
    SUBTRACT_BIAS = "Y"
    SUBTRACT_DARK = "Y"
    DIVIDE_FLAT = "Y"
    CLIP_IMAGE = "Y"
    SUBTRACT_SKY = "Y"
    PLATE_SOLVE = "Y"

    # telescope information
    PIXEL_SIZE = 0.4959  #  0.47  # arcsec per pixel
    NUM_PIXELS = 10560  # pixels per side
    TOROS_DEC_LIMIT = 26.66  # declination limit of the telescope in degrees
    FOV = (PIXEL_SIZE * NUM_PIXELS) / 3600.
    SEARCH_DIST = FOV
    EXP_TIME = 300
    GAIN = 0.380  # in e-/ADU

    # get image information
    AXS_X_RW = 12000
    AXS_Y_RW = 10600
    AXS_X = 10560
    AXS_Y = 10560
    AXS = 10560

    # update the differencing information, primarily the number of stars to use, and the kernel size
    KRNL = 2  # kernel size 2 * KNRL + 1
    STMP = 10  # stamp size ot use 2 * STMP + 1
    ORDR = 1  # order of the kernel to use, 0 is stationary, 1 or 2 is spatially varying
    NRSTARS = 500  # number of stars used to solve for kernel
    BRIGHT_STARS = 20000  # the top stars to search for in kernel stars
    KERNEL_LIMIT = 0.5  # the maximum allowable offset in zeropoint in magnitudes
    AXS_LIMIT = 100  # the number of pixel close to the edge of the frame to use
    RMS_LOW_LIMIT = 0.005  # the lower limit on precision to use for the kernel stars
    RMS_UP_LIMIT = 0.02  # the upper limit on precision to use for the kernel stars

    # update sky subtraction specific information
    PIX = 220

    # a photometry configuration
    FWHM = 15.  # fwhm of the image
    THRESHOLD = 5.  # the threshold for a source above the background

    # aperture information
    APER_SIZE = 16  # circular aperture

    # aperture annulus for the sky background automatically determined from the main aperture
    ANNULI_INNER = APER_SIZE + 2
    ANNULI_OUTER = APER_SIZE + 4

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = "/Users/mcast/Research/toros/ryan/working_directory/"
    ALERTS_DIRECTORY = WORKING_DIRECTORY + 'alerts/'
    ANALYSIS_DIRECTORY = WORKING_DIRECTORY + 'analysis/'
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs/'
    QUERIES_DIRECTORY = WORKING_DIRECTORY + 'queries/'
    CODE_DIFFERENCE_DIRECTORY = WORKING_DIRECTORY + 'difference/'

    # input paths for data etc
    DATA_DIRECTORY = "/Volumes/datadrive/coldharbor/"
    RAW_DIRECTORY = "/Volumes/datadrive/coldharbor/raw/"
    CLEAN_DIRECTORY = DATA_DIRECTORY + "clean/"
    MASTER_MAIN_DIRECTORY = DATA_DIRECTORY + "master/"
    MASTER_DIRECTORY = DATA_DIRECTORY + "master/" + FIELD + "/"
    MASTER_TMP_DIRECTORY = MASTER_DIRECTORY + "tmp_master/"
    CENTROID_DIRECTORY = MASTER_DIRECTORY + "centroids/"
    CALIBRATION_DIRECTORY = DATA_DIRECTORY + "calibration/"
    BIAS_DIRECTORY = CALIBRATION_DIRECTORY + "tmp_bias/"
    FLAT_DIRECTORY = CALIBRATION_DIRECTORY + "tmp_flat/"
    DARK_DIRECTORY = CALIBRATION_DIRECTORY + "tmp_dark/"
    LIGHTCURVE_DIRECTORY = DATA_DIRECTORY + "lc/"
    LIGHTCURVE_FIELD_DIRECTORY = LIGHTCURVE_DIRECTORY + FIELD + "/"
    DIFFERENCED_DIRECTORY = DATA_DIRECTORY + "diff/"
    FLUX_DIRECTORY = DATA_DIRECTORY + "flux/"

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LOG_DIRECTORY, CALIBRATION_DIRECTORY, FLUX_DIRECTORY,
                   QUERIES_DIRECTORY, CLEAN_DIRECTORY, MASTER_MAIN_DIRECTORY, MASTER_DIRECTORY, MASTER_TMP_DIRECTORY,
                   LIGHTCURVE_DIRECTORY, CENTROID_DIRECTORY, RAW_DIRECTORY, BIAS_DIRECTORY, DARK_DIRECTORY,
                   FLAT_DIRECTORY, DIFFERENCED_DIRECTORY, CODE_DIFFERENCE_DIRECTORY, LIGHTCURVE_FIELD_DIRECTORY]

    # BROKER CONFIGURATION SPECIFICS
    LISTEN_NED_WAIT = 1

    # observing conditions
    SEEING = 0.93  # assumes 2 pix FWHM

    # sky brightness at TOLAR in SDSS griz
    SKY = [22.1, 21.1, 20.1, 18.7]

    # SDSS griz bandpass values in nm (width of the filter)
    BP = [147, 141, 147, 147]

    # SDSS griz central wavelength
    CWL = [473.5, 638.5, 775.5, 922.5]

    # telescope information
    TOROS_MIRROR_D = 0.610  # m
    TOROS_MIRROR_R = TOROS_MIRROR_D / 2  # nm

    # CCD information
    READOUT_NOISE = 5.0  # electrons / pixel
    CCD_QE = 0.85
    FILTER_QE = 0.9
    TELESCOPE_SEC_QE = 0.96
    TELECSCOPE_PRI_QE = 0.96
    VIGNETTING = 0.756
    ATMOSPHERE_QE = [0.8, 0.9, 0.9, 0.9]

    # total throughput needs to be multiplied by atmospheric quantum efficiency
    TOTAL_THROUGHPUT = CCD_QE * FILTER_QE * TELESCOPE_SEC_QE * TELECSCOPE_PRI_QE * VIGNETTING

    # force toros field generation?
    FIELD_GENERATION = 'N'

    # The Felix Aguilar Observatory is more Southern than Tolar
    TOROS_LONGITUDE = -69.3265  # -67.32833333
    TOROS_LATITUDE = -31.8023  # -24.62055556
    TOROS_ELEVATION = 2420
    UTC = -3
    MOON_DISTANCE = 60

    EXPOSURE_TIME = 300
    EXPOSURE_TIME_DAY = EXPOSURE_TIME / 60. / 60. / 24.
    NUM_EXPOSURES = 1
    READ_TIME = 90
    OVERHEAD = 30
    TOTAL_EXPOSURE = (EXPOSURE_TIME + READ_TIME) * NUM_EXPOSURES + OVERHEAD
