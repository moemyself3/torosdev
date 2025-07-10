""" This script will prepare images for subtraction and then difference the images from the master frame. """
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from config import Configuration
import os
import numpy as np
from astropy.io import fits
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
from astropy.stats import sigma_clipped_stats
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
import logging
matplotlib.set_loglevel(level = 'warning')
matplotlib.use("TkAgg")
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.INFO)


class BigDiff:

    @staticmethod
    def difference_images():
        """ This function will generate the master frame and generates position files, or if a single frame is chosen,
        then only the position file are generated.

        :return - Nothing is returned, however, the images are differenced
        """

        # get the image list to difference
        files, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY,
                                                     Configuration.FIELD,
                                                     'clean',
                                                     Configuration.FILE_EXTENSION)
        nfiles = len(files)

        # read in the master frame information
        master, master_header = fits.getdata(Configuration.MASTER_DIRECTORY +
                                             Configuration.FIELD +'_master' + Configuration.FILE_EXTENSION, header=True)

        # read in the star list for processing
        star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', delimiter=' ',
                                header=0)

        # prepare the oisdifference.c file for differencing
        BigDiff.prep_ois(master, master_header)

        # begin with the algorithm to difference the images
        for ii in range(0, nfiles):

            fin_nme = Preprocessing.mk_nme(files[ii], 'Y', 'N', 'N', 'N', 'N')

            if os.path.isfile(fin_nme) == 1:
                Utils.log("File " + fin_nme + " found. Skipping...", "info")

            # check to see if the differenced file already exists
            if os.path.isfile(fin_nme) == 0:
                Utils.log("Working to difference file " + files[ii] + ".", "info")
                BigDiff.diff_img(star_list, files[ii], fin_nme)

        Utils.log("Differencing complete for " + Configuration.FIELD + ".", "info")

        return

    @staticmethod
    def diff_img(star_list, file, out_name):
        """ This function will check for and determine the reference stars. It will then difference the image.

        :parameter file - The file name to difference.
        :parameter out_name - The final file name.

        :return - Nothing is returned but the image is differenced
        """

        # read in the image
        org_img, org_header = fits.getdata(file, header=True)

        # read in the master frame header to align the images
        master_header = fits.getheader(Configuration.MASTER_DIRECTORY + Configuration.FIELD + "_master.fits")

        # write the new image file
        img_sbkg = org_img - org_header['SKY']
        img_align = Preprocessing.align_img(img_sbkg, org_header, master_header)

        org_header['WCSAXES'] = master_header['WCSAXES']
        org_header['CRPIX1'] = master_header['CRPIX1']
        org_header['CRPIX2'] = master_header['CRPIX2']
        org_header['PC1_1'] = master_header['PC1_1']
        org_header['PC1_2'] = master_header['PC1_2']
        org_header['PC2_1'] = master_header['PC2_1']
        org_header['PC2_2'] = master_header['PC2_2']
        org_header['CDELT1'] = master_header['CDELT1']
        org_header['CDELT2'] = master_header['CDELT2']
        org_header['CUNIT1'] = master_header['CUNIT1']
        org_header['CUNIT2'] = master_header['CUNIT2']
        org_header['CTYPE1'] = master_header['CTYPE1']
        org_header['CTYPE2'] = master_header['CTYPE2']
        org_header['CRVAL1'] = master_header['CRVAL1']
        org_header['CRVAL2'] = master_header['CRVAL2']
        org_header['LONPOLE'] = master_header['LONPOLE']
        org_header['LATPOLE'] = master_header['LATPOLE']
        org_header['MJDREF'] = master_header['MJDREF']
        org_header['RADESYS'] = master_header['RADESYS']
        org_header['ALIGNED'] = 'Y'

        fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.fits',
                       img_align, org_header, overwrite=True)

        # get the kernel stars for the subtraction
        nstars = BigDiff.find_subtraction_stars_img(img_align, star_list)

        BigDiff.ois_difference(out_name, org_header, nstars)

        return

    @staticmethod
    def ois_difference(out_name, header, nstars):
        """ This function will run the c code oisdifference

        :parameter out_name - The file name for the difference file
        :parameter header - The header of the image
        :parameter nstars - The number of stars in the subtraction

        :return Nothing is returned, the image is differenced
        """
        Utils.log('Now starting image subtraction.', 'info')
        Utils.log('The kernel size is: ' + str(Configuration.KRNL * 2 + 1) + 'x' + str(Configuration.KRNL * 2 + 1) +
                  '; the stamp size is: ' + str(Configuration.STMP * 2 + 1) + 'x' + str(Configuration.STMP * 2 + 1) +
                  '; the polynomial is: ' + str(Configuration.ORDR) + ' order; and ' + str(nstars) +
                  ' stars were used in the subtraction.', 'info')

        # change to the directory
        os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)

        # run the c code
        shh = os.system('./a.out')

        # update the header file
        dimg, diff_header = fits.getdata('dimg.fits', header=True)
        header['diffed'] = 'Y'
        header['nstars'] = nstars

        # now mask the missing master frame parts #### THIS WILL CHANGE PER FIELD!!!! LIKELY YOU SHOULD REMOVE#####
        dimg[0:490, :] = 0
        dimg[10045:-1, :] = 0
        dimg[:, 0:530] = 0
        dimg[:, 10465:-1] = 0
        #### END LIKELY REMOVE

        # update the image with the new file header
        fits.writeto('dimg.fits', dimg, header, overwrite=True)

        # move the differenced file to the difference directory
        os.system('mv dimg.fits ' + out_name)

        # change back to the working directory
        os.chdir(Configuration.WORKING_DIRECTORY)

        # get the photometry from the differenced image
        Utils.log('Image subtraction complete.', 'info')
        return

    @staticmethod
    def prep_ois(master, master_header):
        """ This function will prepare the files necessary for the ois difference.

        :parameter master - The master image for differencing
        :parameter master_header - The header file for the master image

        :return - Nothing is returned but the necessary text files are written,
                    and the code is compiled for differencing
        """
        # compile the oisdifference.c code
        os.system('cp oisdifference.c ' + Configuration.CODE_DIFFERENCE_DIRECTORY)
        os.chdir(Configuration.CODE_DIFFERENCE_DIRECTORY)
        os.system('gcc `pkg-config --cflags --libs cfitsio` oisdifference.c')
        os.chdir(Configuration.WORKING_DIRECTORY)

        # write the new master file
        fits.writeto(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.fits', master, master_header, overwrite=True)

        # prepare the text files
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'ref.txt', 'w', "ref.fits")
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'img.txt', 'w', "img.fits")

        return

    @staticmethod
    def find_subtraction_stars_img(img, star_list):
        """ This function will find the subtraction stars to use for the differencing, they will be the same stars for
        every frame. This will help in detrending later.

        """

        Utils.log('Finding stars for kernel from the star list.', 'info')
        diff_list = star_list.copy().reset_index(drop=True)

        ## REMOVE THIS PART FOR NON-47TUC FIELDS!!!!
        # remove stars near 47 Tuc and the small cluster

        diff_list['bd_star'] = np.where((diff_list['xcen'] > 4300) & (diff_list['xcen'] < 9300) &
                                        (diff_list['ycen'] > 3600) & (diff_list['ycen'] < 8200), 1, 0)
        diff_list['bd_star'] = np.where((diff_list['xcen'] > 1200) & (diff_list['xcen'] < 1900) &
                                         (diff_list['ycen'] > 5100) & (diff_list['ycen'] < 5600), 2,
                                        diff_list['bd_star'])
        diff_list = diff_list[diff_list.bd_star == 0].copy().reset_index(drop=True)
        ## END REMOVAL FOR NON-47TUC FIELDS!!!!

        # now check for stars based on their magnitude differences
        positions = np.transpose((diff_list['xcen'], diff_list['ycen']))
        aperture = CircularAperture(positions, r=Configuration.APER_SIZE)

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, aperture, method='exact')

        # the global background was subtracted, and teh master mangitude does not have exposure time corrected
        flux = np.array(phot_table['aperture_sum']) * Configuration.GAIN

        # convert to magnitude
        mag = 25 - 2.5 * np.log10(flux)

        # clip likely variables, these objects have large magnitudes changes relative to the master frame
        diff_list['dmag'] = diff_list['master_mag'].to_numpy() - mag

        # the clip will be 2 sigma-clipping above and below the mean offset
        dmn, dmd, dsg = sigma_clipped_stats(diff_list.dmag, sigma=2)
        dmag_plus = dmd + dsg
        dmag_minus = dmd - dsg
        diff_list = diff_list[(diff_list['dmag'] < dmag_plus) &
                                 (diff_list['dmag'] > dmag_minus)].copy().reset_index(drop=True)

        # now check for non-crowded stars
        diff_list['prox'] = diff_list.apply(lambda x: np.sort(np.sqrt((x.xcen - diff_list.xcen) ** 2 +
                                                                          (x.ycen - diff_list.ycen) ** 2))[1], axis=1)
        diff_list = diff_list[diff_list.prox > (2 * Configuration.STMP + 1)].copy().reset_index(drop=True)

        # make a magnitude cut
        diff_list = diff_list[(diff_list.phot_g_mean_mag < 15) &
                              (diff_list.master_mag > 10)].copy().reset_index(drop=True)

        if len(diff_list) > Configuration.NRSTARS:
            diff_list = diff_list.sample(n=Configuration.NRSTARS)
            nstars = Configuration.NRSTARS
        else:
            nstars = len(diff_list)
            Utils.log("There are not enough stars on the frame to use " + str(Configuration.NRSTARS) +
                      " in the subtraction. Instead, all available stars (" +
                      str(nstars) + ") will be used in the subtraction.",
                      "info")

        # add 1 for indexing in C vs indexing in python
        diff_list['x'] = np.around(diff_list['xcen'] + 1, decimals=0)
        diff_list['y'] = np.around(diff_list['ycen'] + 1, decimals=0)

        # write the parameter file for the C code to use
        Utils.write_txt(Configuration.CODE_DIFFERENCE_DIRECTORY + 'parms.txt', 'w', "%1d %1d %1d %4d\n" %
                        (Configuration.STMP, Configuration.KRNL, Configuration.ORDR, nstars))

        # export the differencing stars
        diff_list[['x', 'y']].astype(int).to_csv(Configuration.CODE_DIFFERENCE_DIRECTORY +
                                                 'refstars.txt', index=0, header=0, sep=" ")

        return nstars
