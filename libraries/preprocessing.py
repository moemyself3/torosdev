import pandas as pd
from config import Configuration
from libraries.utils import Utils
import os
import numpy as np
import astropy.wcs as pywcs
import scipy.ndimage
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import twirl
from photutils.detection import DAOStarFinder
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats
#import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt

class Preprocessing:

    @staticmethod
    def mk_combined_bias_and_dark(image_overwrite):
        """ This function will make a master bias frame (sans a mean pixel value) and a master dark frame. The master
        dark frame uses the bias frames from the same night to generate. The master bias frame has its median pixel
        value removed so it can be scaled to the value of the overscan region later. This is because the TOROS bias
        level seems to be changing with time.

        :parameter image_overwrite - Y/N if you want to force an overwrite for each master frame

        :return - The master bias and the master dark are returned
        """

        # set the names of the master frames
        bias_name = 'bias' + Configuration.FILE_EXTENSION
        dark_name = 'dark' + Configuration.FILE_EXTENSION

        # check to see if the master bias frame exists, and if it does save some time by skipping
        if ((os.path.isfile(Configuration.CALIBRATION_DIRECTORY + dark_name) == 0) |
                (os.path.isfile(Configuration.CALIBRATION_DIRECTORY + bias_name) == 0) |
                (image_overwrite == 'Y')):

            # get the image lists
            biases = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'bias_list.csv', sep=',')
            darks = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'dark_list.csv', sep=',')
            total_bias = 0
            total_dark = 0
            # group based on the date to get the unique dates
            bias_dates = biases.groupby('Date').agg({'count'}).reset_index()['Date'].to_list()
            ndates = len(bias_dates)

            # make holders for the total number of nights
            bias_bulk = np.ndarray(shape=(ndates, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))
            dark_bulk = np.ndarray(shape=(ndates, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))
            bias_bulk_filepath = []
            dark_bulk_filepath = []
            zdx = 0
            for dte in bias_dates:
                # determine how many bias and dark frames exist on this date
                bias_frames = biases[biases.Date == dte]
                bias_list = bias_frames.apply(lambda x: Configuration.DATA_DIRECTORY + "bias/"  + x.Date + x.File, axis=1).to_list()
                nbias = len(bias_frames)

                dark_frames = darks[darks.Date == dte]
                ndarks = len(dark_frames)

                # if dark frames exist, then move forward
                if ndarks > 0:
                    Utils.log("Working to make mini bias and dark files for " + dte + '.', "info")
                    dark_list = dark_frames.apply(lambda x: Configuration.DATA_DIRECTORY + "darks/" + x.Date + x.Files,
                                                  axis=1).to_list()

                    # generate the frame holder
                    bias_hld = np.ndarray(shape=(nbias, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))
                    dark_hld = np.ndarray(shape=(ndarks, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))

                    idx = 0
                    for ii in range(0, nbias):
                        # read in the bias file
                        bias_tmp, bias_head = fits.getdata(bias_list[ii], header=True)
                        bias_hld[idx] = bias_tmp

                        idx = idx + 1
                        total_bias = total_bias + 1
                        bias_tmp = None
                    bias_tmp_mdn = np.median(bias_hld, axis=0)
                    bias_hld = None
                    Utils.log("Bias done.", "info")

                    jdx = 0
                    for jj in range(0, ndarks):
                        # read in the dark file
                        dark_tmp, dark_head = fits.getdata(dark_list[jj], header=True)

                        # remove the bias frame from the temporary dark
                        dark_hld[jdx] = dark_tmp - bias_tmp_mdn
                        jdx = jdx + 1
                        total_dark = total_dark + 1
                        dark_tmp = None
                    dark_tmp_mdn = np.median(dark_hld, axis=0)
                    dark_hld = None
                    # write the image out to the temporary directory
                    dark_tmp_filename = Configuration.DARK_DIRECTORY + str(zdx).zfill(2) + "_scale_tmp_dark.fits"
                    fits.writeto(dark_tmp_filename,
                                 dark_tmp_mdn, overwrite=True)
                    Utils.log("Dark done.", "info")

                    # what are the sizes of the over scan?
                    ovs_x_size = 180
                    ovs_y_size = 20  # it is 20 above and 20 below

                    # what are the sizes of each chip (not including the over scan)?
                    chip_x_size = 1320
                    chip_y_size = 5280

                    # what is the full size of the chip (including over scan)
                    full_chip_x = chip_x_size + ovs_x_size
                    full_chip_y = chip_y_size + ovs_y_size

                    # move through x and y to mask the "image" parts of hte image
                    for x in range(0, Configuration.AXS_X_RW, full_chip_x):
                        for y in range(0, Configuration.AXS_Y_RW, full_chip_y):

                            # put the clipped image into the holder image
                            bias_tmp_mdn[y:y + full_chip_y, x:x + full_chip_x] = (
                                    bias_tmp_mdn[y:y + full_chip_y, x:x + full_chip_x] -
                                    np.median(bias_tmp_mdn[y:y + full_chip_y, x:x + full_chip_x]))
                    bias_tmp_filename = Configuration.BIAS_DIRECTORY + str(zdx).zfill(2) + "_scale_tmp_bias.fits"
                    fits.writeto(bias_tmp_filename,
                                 bias_tmp_mdn, overwrite=True)

                    # update the bulk holder
                    #bias_bulk[zdx] = bias_tmp_mdn
                    #dark_bulk[zdx] = dark_tmp_mdn
                    bias_bulk_filepath.append(bias_tmp_filename)
                    dark_bulk_filepath.append(dark_tmp_filename)
                    zdx = zdx + 1
                    bias_tmp_mdn = None
                    dark_tmp_mdn = None

            # load bias bulk
            for index, tmp_bias_file in enumerate(bias_bulk_filepath):
                bias_bulk[index] = fits.getdata(tmp_bias_file)

            # update the header with relevant information
            bias_hdu = fits.PrimaryHDU()
            bias_header = bias_hdu.header
            bias_header['BIAS_COMB'] = 'median'
            bias_header['NUM_BIAS'] = total_bias
            bias = np.median(bias_bulk[0:zdx], axis=0)
            bias_bulk = None
            # write the image out to the master directory
            fits.writeto(Configuration.CALIBRATION_DIRECTORY + bias_name,
                         bias, bias_header, overwrite=True)

            # load dark bulk
            for index, tmp_dark_file in enumerate(dark_bulk_filepath):
                dark_bulk[index] = fits.getdata(tmp_dark_file)

            # update the header with relevant information
            dark_hdu = fits.PrimaryHDU()
            dark_header = dark_hdu.header
            dark_header['DARK_COMB'] = 'median'
            dark_header['NUM_DARK'] = total_dark
            dark_header['EXPTIME'] = 300
            dark_header['BIAS_SUB'] = 'Y'
            dark = np.median(dark_bulk[0:zdx], axis=0)
            dark_bulk = None
            # write the image out to the master directory
            fits.writeto(Configuration.CALIBRATION_DIRECTORY + dark_name,
                         dark, dark_header, overwrite=True)

            Utils.log("Scalable bias and dark frames generated.", "info")

        else:
            bias = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'bias.fits')
            dark = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'dark.fits')
        return bias, dark

    @staticmethod
    def mk_flat(flat_exp=5, dark_exp=300):
        """ This function will make the master flat frame using the provided image list.
        :parameter flat_exp - The exposure time for the flat frames
        :parameter dark_exp - The exposure time for the dark frames

        :return - The flat field for the given date is returned and written to the calibration directory
        """

        if os.path.isfile(Configuration.CALIBRATION_DIRECTORY + "flat.fits") == 0:

            # read in the bias frame and dark frame
            bias = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'bias.fits')
            dark = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'dark.fits')

            # what are the sizes of the over scan?
            ovs_x_size = 180
            ovs_y_size = 20

            # what are the sizes of each chip (not including the over scan)?
            chip_x_size = 1320
            chip_y_size = 5280

            # what is the full size of the chip (including over scan)
            full_chip_x = chip_x_size + ovs_x_size
            full_chip_y = chip_y_size + ovs_y_size

            # get the image list
            images = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'flat_list.csv', sep=',')
            image_list = images.apply(lambda x: Configuration.DATA_DIRECTORY + "flats/" + x.Date + x.File, axis=1).to_list()

            # determine the number of loops we need to move through for each image
            nfiles = len(image_list)
            nbulk = 20

            # get the integer and remainder for the combination
            full_bulk = nfiles // nbulk
            part_bulk = nfiles % nbulk

            if part_bulk > 0:
                hold_bulk = full_bulk + 1
            else:
                hold_bulk = full_bulk

            # here is the 'holder'
            hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))
            hold_data_filepath = []

            # update the log
            Utils.log("Generating a master flat field from multiple files in bulks of " + str(nbulk) +
                      " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                      str(hold_bulk) + " mini-files to combine.", "info")

            tmp_num = 0
            for kk in range(0, hold_bulk):

                # loop through the images in sets of nbulk
                if kk < full_bulk:
                    # generate the image holder
                    block_hold = np.ndarray(shape=(nbulk, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))

                    # generate the max index
                    mx_index = nbulk
                else:
                    # generate the image holder
                    block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS_Y_RW, Configuration.AXS_X_RW))

                    # generate the max index
                    mx_index = part_bulk

                # make the starting index
                loop_start = kk * nbulk
                idx_cnt = 0

                Utils.log("Making mini-flat field frame number " + str(kk) + ".", "info")

                # now loop through the images
                for jj in range(loop_start, mx_index + loop_start):
                    # make a copy of the bias frame
                    bias_scl = bias.copy()

                    # read in the flat file
                    flat_tmp, flat_head = fits.getdata(image_list[jj], header=True)

                    # get the scale factor for the dark frame
                    dark_scale = dark_exp / flat_exp

                    # move through x and y to mask the "image" parts of hte image
                    for x in range(0, Configuration.AXS_X_RW, full_chip_x):
                        for y in range(0, Configuration.AXS_Y_RW, full_chip_y):

                            if y == 0:
                                # pull out the overscan from the raw image and set the "image" to 0
                                img_slice = flat_tmp[y:y + full_chip_y, x:x + full_chip_x].copy()
                                img_slice[0:chip_y_size, 0:chip_x_size] = 0
                            else:
                                # pull out the overscan from the raw image and set the "image" to 0
                                img_slice = flat_tmp[y:y + full_chip_y, x:x + full_chip_x].copy()
                                img_slice[ovs_y_size:ovs_y_size + chip_y_size, 0:chip_x_size] = 0

                            # put the clipped image into the holder image
                            bias_scl[y:y + full_chip_y, x:x + full_chip_x] = (
                                        bias_scl[y:y + full_chip_y, x:x + full_chip_x] +
                                        np.median(img_slice[img_slice > 0]))
                            img_slice = None

                    # now remove the bias from the image
                    flat_bias = flat_tmp - bias_scl

                    # now remove the dark current from the image
                    block_hold[idx_cnt] = flat_bias - (dark / dark_scale)

                    if (jj % 20) == 0:
                        Utils.log(str(jj + 1) + " files read in. " +
                                  str(nfiles - jj - 1) + " files remain.", "info")

                    # increase the iteration
                    idx_cnt += 1
                    flat_tmp = None # set to None since we use it again
                    flat_bias = None
                    bias_scl = None
                # median the data into a single file
                #hold_data[kk] = np.median(block_hold, axis=0)
                block_hold_median = np.median(block_hold, axis=0)
                block_hold = None

                # write out the temporary file
                flat_tmp_filename = Configuration.FLAT_DIRECTORY + str(tmp_num).zfill(2) + "_tmp_flat.fits"
                fits.writeto(flat_tmp_filename,
                             block_hold_median, overwrite=True)
                hold_data_filepath.append(flat_tmp_filename)
                block_hold_median = None
                tmp_num = tmp_num + 1

            # load tmp files into hold_data
            for index, tmp_file in enumerate(hold_data_filepath):
                hold_data[index] = fits.getdata(tmp_file)

            # median the mini-images into one large image
            flat_image = np.median(hold_data, axis=0)
            nflat_image = flat_image / np.median(flat_image[5950:5950 + 3900, 3200:3200 + 900])

            # pull the header information from the first file of the set
            flat_header = fits.getheader(image_list[0])
            flat_header['comb_typ'] = 'median'
            flat_header['median_val'] = np.median(flat_image[5950:5950 + 3900, 3200:3200 + 900])
            flat_header['norm_pix'] = 'median'
            flat_header['num_comb'] = len(image_list)
            flat_header['mean_pix'] = np.mean(nflat_image)
            flat_header['std_pix'] = np.std(nflat_image)
            flat_header['max_pix'] = np.max(nflat_image)
            flat_header['min_pix'] = np.min(nflat_image)
            flat_header['mean'] = np.mean(flat_image)
            flat_header['std'] = np.std(flat_image)
            flat_header['max'] = np.max(flat_image)
            flat_header['min'] = np.min(flat_image)

            # write the image out to the master directory
            fits.writeto(Configuration.CALIBRATION_DIRECTORY + "flat.fits",
                         nflat_image, flat_header, overwrite=True)

        else:
            nflat_image = fits.getdata(Configuration.CALIBRATION_DIRECTORY + "flat.fits")

        return nflat_image

    @staticmethod
    def correct_header(img, header):
        """ This function will plate solve the image, add the time stamp, and exposure time to the header if need be.

        :parameter img - This is the image you want to plate solve
        :parameter header - the header of the image you want to correct

        return img, header - The corrected image will be sent back
        """

        # get the approximate center of the image
        center = SkyCoord(Configuration.RA, Configuration.DEC, unit=['deg', 'deg'])
        pixel = Configuration.PIXEL_SIZE * u.arcsec
        fov = np.max(np.shape(img)) * pixel.to(u.deg)

        # now query the gaia region
        all_stars = twirl.gaia_radecs(center, 1.25 * fov)

        # keep the isolated stars
        all_stars = twirl.geometry.sparsify(all_stars, 0.01)[0:30]

        # get the stars in the image
        mean, median, std = sigma_clipped_stats(img, sigma=3.0)
        daofind = DAOStarFinder(fwhm=3.0, threshold=100)
        sources = daofind(img - median)

        # sort based on "real" stars and only select the top 50 or so
        sources = sources[(sources['flux'] > 0) &
                          (sources['xcentroid'] > 2500) &
                          (sources['ycentroid'] > 2500)]
        sources.sort('flux', reverse=True)
        sources = sources[0:100]

        # generate the set up for twirl to match stars
        xy = np.ndarray(shape=(len(sources), 2))
        idx = 0
        for xx, yy in zip(sources['xcentroid'], sources['ycentroid']):
            xy[idx][0] = xx
            xy[idx][1] = yy
            idx = idx + 1

        # remove double stars
        xy = twirl.geometry.sparsify(xy, 30)

        # now begin to match stars and make the header
        try:
            # now compute the new wcs
            wcs = twirl.compute_wcs(xy, all_stars)

            # add the WCS to the header
            h = wcs.to_header()

            # transfer the header information
            for idx, v in enumerate(h):
                 header[v] = (h[idx], h.comments[idx])

            # in some cases twirl will fail, but still complete, mark this in the image header so we are aware
            if header['PC1_1'] < -1:
                header['BAD_WCS'] = 'Y'
            else:
                header['BAD_WCS'] = 'N'

        except:
            Utils.log("Bad image!", "info")
            header['BADIMAGE'] = 'Y'

        return img, header

    @staticmethod
    def sky_subtract(img, header, sky_write='N'):
        """  The function has been updated to include the photutils background subtraction routine.

        :parameter img - The image to be cleaned
        :parameter header - The header object to be updated
        :parameter sky_write - Y/N if you want to write the residual background for de-bugging

        :return img_sub, header - The cleaned image and updated header file
        """
        # set up the background clipping functions
        sigma_clip = SigmaClip(sigma=3)
        bkg_estimator = MedianBackground()

        #### REMOVE THIS PART FOR NON-47TUC FIELDS!!!!
 #       # find the most likely position of the cluster
 #       daofind = DAOStarFinder(fwhm=3.0, threshold=50)
 #       sources = daofind(img[3000:, 3000:])

 #       # find the cluster
 #       x_cen = (np.sum(sources[sources['flux'] > 0]['xcentroid'] * sources[sources['flux'] > 0]['flux']) /
 #                np.sum(sources[sources['flux'] > 0]['flux'])) + 3000
 #       y_cen = (np.sum(sources[sources['flux'] > 0]['ycentroid'] * sources[sources['flux'] > 0]['flux']) /
 #                np.sum(sources[sources['flux'] > 0]['flux'])) + 3000

 #       # make the masked image
 #       mask_img = np.zeros((Configuration.AXS_Y, Configuration.AXS_X))

 #       mask_img[int(y_cen - 1000):int(y_cen + 1000), int(x_cen - 1000):int(x_cen + 1000)] = 1  # TUC-47
 #       if ((x_cen - 5200) - 200 > 0) & ((x_cen - 5200) > 0):
 #           mask_img[int((y_cen - 500) - 200):int((y_cen - 500) + 200),
 #           int((x_cen - 5200) - 200):int((x_cen - 5200) + 200)] = 1 # Other GC
 #       #### END REMOVAL FOR NON-47TUC FIELDS!!!!
        mask_img = None
        # do the 2D background estimation, if there is no mask, then remove mask_img
        bkg = Background2D(img, (Configuration.PIX, Configuration.PIX), filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask_img)
        sky = bkg.background

        # subtract the sky gradient and add back the median background
        img_sub = img - sky
        fin_img = img_sub + bkg.background_median

        # now correct the pixel line near the overscan
        x_skip = 1320
        y_skip = 220
        y_full = 5280
        for x in range(0, Configuration.AXS_X, x_skip):
            for y in range(0, Configuration.AXS_Y, y_skip):

                # get the sky from the slice and full y column
                mdn_slc = np.median(fin_img[y:y + y_skip, x])
                mdn = np.median(fin_img[y:y + y_full, x])

                # update the column with the new background, unless something chonky is thowing off the calculation
                # like 47-Tuc or a galaxy etc
                if mdn_slc < (mdn + 0.1 * mdn):
                    fin_img[y:y + y_skip, x] = (fin_img[y:y + y_skip, x] -
                                                np.median(fin_img[y:y + y_skip, x]) +
                                                bkg.background_median)
                else:
                    fin_img[y:y + y_skip, x] = (fin_img[y:y + y_skip, x] -
                                                np.median(fin_img[y:y + y_full, x]) +
                                                bkg.background_median)

        # update the header
        header['sky'] = bkg.background_median
        header['sky_sig'] = bkg.background_rms_median
        header['sky_sub'] = 'yes'

        # if desired, write out the sky background to the working directory
        if sky_write == 'Y':
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'sky_background' + Configuration.FILE_EXTENSION,
                         sky, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img' + Configuration.FILE_EXTENSION,
                         img, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img_sub' + Configuration.FILE_EXTENSION,
                         fin_img, header=header, overwrite=True)

        return fin_img, header

    @staticmethod
    def subtract_scaled_bias_dark(img, header):
        """ This function will scale the bias frame to match the overscan and then remove the dark level.

        :parameter - img - The image to remove the bias and dark frame frome
        :parameter - header - The header of the image

        :return img_bias_dark, header
        """

        # read in the bias frame and dark frame
        bias = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'bias.fits')
        dark = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'dark.fits')

        # what are the sizes of the over scan?
        ovs_x_size = 180
        ovs_y_size = 20

        # what are the sizes of each chip (not including the over scan)?
        chip_x_size = 1320
        chip_y_size = 5280

        # what is the full size of the chip (including over scan)
        full_chip_x = chip_x_size + ovs_x_size
        full_chip_y = chip_y_size + ovs_y_size

        # make a copy of the bias frame
        bias_scl = bias.copy()
        del bias # del bias because we do not use it again

        # move through x and y to mask the "image" parts of hte image
        for x in range(0, Configuration.AXS_X_RW, full_chip_x):
            for y in range(0, Configuration.AXS_Y_RW, full_chip_y):

                if y == 0:
                    # pull out the overscan from the raw image and set the "image" to 0
                    img_slice = img[y:y + full_chip_y, x:x + full_chip_x].copy()
                    img_slice[0:chip_y_size, 0:chip_x_size] = 0
                else:
                    # pull out the overscan from the raw image and set the "image" to 0
                    img_slice = img[y:y + full_chip_y, x:x + full_chip_x].copy()
                    img_slice[ovs_y_size:ovs_y_size + chip_y_size, 0:chip_x_size] = 0

                # put the clipped image into the holder image
                bias_scl[y:y + full_chip_y, x:x + full_chip_x] = (bias_scl[y:y + full_chip_y, x:x + full_chip_x] +
                                                                  np.median(img_slice[img_slice > 0]))

        # now remove the bias from the image
        img_bias = img - bias_scl

        # now remove the dark current from the image
        img_bias_dark = img_bias - dark

        # now update the image header
        # update the header
        header['BIAS_SUBT'] = 'Y'
        header['BIAS_TYPE'] = 'SCALE'
        header['DARK_SUBT'] = 'Y'

        return img_bias_dark, header

    @staticmethod
    def clip_image(img, header):
        """ This function will clip the image and remove any overscan regions. This function is written for TOROS
        specifically, and you will need to update it for any given CCD.

        :parameter - image - The image to clip
        :parameter - header - the header of the image

        :return image_clip, header - The clipped image and the new header
        """

        # make the clipped image
        image_clip = np.zeros((Configuration.AXS_Y, Configuration.AXS_X))

        # what are the sizes of the over scan?
        ovs_x_size = 180
        ovs_y_size = 20

        # what are the sizes of each chip (not including the over scan)?
        chip_x_size = 1320
        chip_y_size = 5280

        # what is the full size of the chip (including over scan)
        full_chip_x = chip_x_size + ovs_x_size
        full_chip_y = chip_y_size + ovs_y_size

        # move through x and y
        idx = 0
        for x in range(0, Configuration.AXS_X_RW, full_chip_x):

            idy = 0
            for y in range(0, Configuration.AXS_Y_RW, full_chip_y):

                if y == 0:
                    # put the clipped image into the holder image
                    image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y:y + chip_y_size, x:x + chip_x_size]
                else:
                    # put the clipped image into the holder image
                    image_clip[idy:idy + chip_y_size, idx:idx + chip_x_size] = img[y + ovs_y_size:y + ovs_y_size + chip_y_size, x:x + chip_x_size]

                # increase the size of the yclip
                idy = idy + chip_y_size

            # increase the size of the xclip
            idx = idx + chip_x_size

        # update the header
        header['OVERSCAN'] = 'removed'
        header['X_CLIP'] = ovs_x_size
        header['Y_CLIP'] = ovs_y_size

        return image_clip, header

    @staticmethod
    def flat_divide(img, header):
        """ This function will divide a flat field.
        :parameter img - The image to flatten
        :parameter header - The image header file

        :return flat_div, header - The updated image and header """
        # read in the flat frame
        flat = fits.getdata(Configuration.CALIBRATION_DIRECTORY + "flat.fits")

        # subtract the bias from the image
        flat_div = img / flat

        # update the header
        header['FLATTEN'] = 'Y'

        return flat_div, header

    @staticmethod
    def mk_nme(file, difference_image='N', image_clip='N', sky_subtract='N', bias_subtract='N',
               flat_divide='N',  dark_subtract="N", plate_solve='N'):
        """ This function will create the appropriate name for the file based on while steps are taken.
        :parameter file - The string with the file name
        :parameter image_clip - Y/N if the image was clipped
        :parameter bias_subtract - Y/N if a bias is subtracted
        :parameter flat_divide - Y/N if a flat field is divided
        :parameter dark_subtract - Y/N if the image was dark subtracted
        :parameter sky_subtract - Y/N if sky subtraction was taken
        :parameter difference_image - Y/N if image subtraction occurred
        :parameter plate_solve - Y/N if the plate solving occurred

        :return file_name - A string with the new file name
        """

        # if everything is N then the file name is the original filename
        file_name = file

        # update the file name with a 'd' if at the differencing step
        if difference_image == 'Y':
            file = file_name.replace("/clean/", "/diff/")
            nme_hld = file.split('.fits')
            file_name = nme_hld[0]  + 'ad' + Configuration.FILE_EXTENSION

        # otherwise...
        if difference_image == 'N':
            # first replace the "raw" directory with the "clean" directory
            file_hld = file_name.split('/')
            file = Configuration.CLEAN_DIRECTORY + file_hld[-2] + "/" + Configuration.FIELD + "/" + file_hld[-1]
            nme_hld = file.split('.fits')

            # update the name to be appropriate for what was done to the file
            # nothing occurs
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'N') \
                    and (dark_subtract == 'N') and (image_clip == 'N'):
                file_name =  nme_hld[0]  + Configuration.FILE_EXTENSION
            # bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_b' + Configuration.FILE_EXTENSION
            # flat
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_f' + Configuration.FILE_EXTENSION
            # sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_s' + Configuration.FILE_EXTENSION
            # dark subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_k' + Configuration.FILE_EXTENSION
            # image clipping only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_c' + Configuration.FILE_EXTENSION
            
            # bias and clip only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bc' + Configuration.FILE_EXTENSION
            # bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bf' + Configuration.FILE_EXTENSION
            # bias and sky_subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bs' + Configuration.FILE_EXTENSION
            # bias and dark only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bk' + Configuration.FILE_EXTENSION
                
            # bias and flat and sky subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bfs' + Configuration.FILE_EXTENSION
            # bias and flat and dark
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bkf' + Configuration.FILE_EXTENSION
            # bias and flat and image clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bcf' + Configuration.FILE_EXTENSION

            # bias and flat and sky and dark
            if (bias_subtract == 'Y') and (flat_divide == 'Y')and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='N'):
                file_name =  nme_hld[0]  + '_bkfs' + Configuration.FILE_EXTENSION
            # bias and flat and sky and clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y')and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bcfs' + Configuration.FILE_EXTENSION

            # bias and flat and dark and sky and clip
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='Y'):
                file_name =  nme_hld[0]  + '_bkcfs' + Configuration.FILE_EXTENSION

            # bias and flat and dark and sky and clip and plate_solve
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y') and (image_clip =='Y') and (plate_solve == 'Y'):
                file_name =  nme_hld[0]  + '_bkfcsp' + Configuration.FILE_EXTENSION
                
        return file_name

    @staticmethod
    def align_img(image, header1, header2, preserve_bad_pixels=False):
        """
        This function is based on the FITS_tools utility hcongrid and it interpolates an image from one .fits header
        to another

        :parameter image - the image to transform
        :parameter header1 - the header of the image
        :parameter header2 - the header to transform to
        :parameter preserve_bad_pixels - Try to set NAN pixels to NAN in the zoomed image.  Otherwise, bad
            pixels will be set to zero

        :return new_image - The image which has been transformed to the reference header

        """

        # convert the headers to WCS objects
        wcs1 = pywcs.WCS(header1)
        wcs1.naxis1 = header1['NAXIS1']
        wcs1.naxis2 = header1['NAXIS2']
        wcs2 = pywcs.WCS(header2)
        wcs2.naxis1 = header2['NAXIS1']
        wcs2.naxis2 = header2['NAXIS2']

        # get the shape
        outshape = [wcs2.naxis2, wcs2.naxis1]
        yy2, xx2 = np.indices(outshape)

        # get the world coordinates of the output image
        lon2, lat2 = wcs2.wcs_pix2world(xx2, yy2, 0)
        xx1, yy1 = wcs1.wcs_world2pix(lon2, lat2, 0)

        # make a grid for the image to transform to
        grid1 = np.array([yy1.reshape(outshape), xx1.reshape(outshape)])

        # identify bad pixels
        bad_pixels = np.isnan(image) + np.isinf(image)
        image[bad_pixels] = 0

        # make the new image
        newimage = scipy.ndimage.map_coordinates(image, grid1)

        # replace bad pixels
        if preserve_bad_pixels:
            newbad = scipy.ndimage.map_coordinates(bad_pixels, grid1, order=0,
                                                   mode='constant',
                                                   cval=np.nan)
            newimage[newbad] = np.nan

        return newimage
