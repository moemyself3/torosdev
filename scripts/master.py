from config import Configuration
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
import numpy as np
import os
from astropy.io import fits
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.centroids import centroid_sources
import pandas as pd
import warnings
from astroquery.mast import Catalogs
from astropy.wcs import WCS
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats


class Master:

    @staticmethod
    def pull_master():
        """ This script will generate the master file and photometry file for the data reduction.

        :return - The master frame is returned and the star list is printed
        """

        master, master_header = Master.mk_master()

        star_list = Master.master_phot(master, master_header)

        return master, star_list  # kernel_stars

    @staticmethod
    def master_phot(master, master_header):
        """ This program will generate the star lists for the master frame and provide a photometry file."""

        if os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt') == 0:

            # create the string useful for query_region
            field = str(Configuration.RA) + " " + str(Configuration.DEC)

            # select the columns we want to import into the data table
            columns = ["toros_field_id", "source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag",
                       "phot_rp_mean_mag", "teff_val", "parallax", "parallax_error", "pmra",
                       "pmra_error", "pmdec", "pmdec_error"]

            # run the query
            Utils.log('Querying MAST for all stars within the toros field: ' + str(Configuration.FIELD), 'info')
            catalog_data = Catalogs.query_region(field,
                                                 radius=Configuration.SEARCH_DIST,
                                                 catalog="Gaia").to_pandas()
            Utils.log('Query finished. ' + str(len(catalog_data)) + ' stars found.', 'info')

            # add the toros field to the catalog data
            catalog_data['toros_field_id'] = Configuration.FIELD

            # pull out the necessary columns
            star_list = catalog_data[columns]

            # get the header file and convert to x/y pixel positions
            w = WCS(master_header)
            ra = star_list.ra.to_numpy()
            dec = star_list.dec.to_numpy()

            # convert to x, y
            x, y = w.all_world2pix(ra, dec, 0)

            # add the x/y to the star data frame
            star_list['x'] = x
            star_list['y'] = y
            star_list = star_list[(star_list.x >= 530) & (star_list.x < 10465) &
                                  (star_list.y >= 490) & (star_list.y < 10045)].copy().reset_index(drop=True)

            star_list['xcen'], star_list['ycen'] = centroid_sources(master,
                                                                    star_list.x.to_numpy(),
                                                                    star_list.y.to_numpy(),
                                                                    box_size=5)
            bd_idxs = np.where(np.isnan(star_list.xcen) | np.isnan(star_list.ycen))
            if len(bd_idxs[0]) > 0:
                for bd_idx in bd_idxs[0]:
                    star_list.loc[bd_idx, 'xcen'] = star_list.loc[bd_idx, 'x']
                    star_list.loc[bd_idx, 'ycen'] = star_list.loc[bd_idx, 'y']

            # centroid the positions
            positions = star_list[['xcen', 'ycen']].copy().reset_index(drop=True)  # positions = (x, y)

            # run aperture photometry
            # set up the star aperture and sky annuli
            aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
            annulus_aperture = CircularAnnulus(positions, r_in=Configuration.ANNULI_INNER,
                                               r_out=Configuration.ANNULI_OUTER)
            apers = [aperture, annulus_aperture]

            # run the photometry to get the data table
            phot_table = aperture_photometry(master, apers, method='exact')

            # extract the sky background for each annuli based on either a global or local subtraction
            sky = phot_table['aperture_sum_1'] / annulus_aperture.area

            # subtract the sky background to get the stellar flux and square root of total flux to get the error
            flux = np.array(phot_table['aperture_sum_0'])

            # calculate the expected photometric error
            flux_er = np.sqrt((phot_table['aperture_sum_0']))

            # convert to magnitude
            mag = 25. - 2.5 * np.log10(flux)
            mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

            # initialize the light curve data frame
            star_list['master_mag'] = mag
            star_list['master_mag_er'] = mag_er
            star_list['master_flux'] = flux
            star_list['master_flux_er'] = flux_er
            star_list['sky'] = sky
            star_list = star_list[star_list['master_flux'] > 0]
            star_list = star_list.sort_values(by='master_mag').reset_index()
            star_list = star_list.rename(columns={'index': 'star_id'})

            star_list.to_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', sep=' ', index=False)

        else:
            star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt', delimiter=' ', header=0)

        return star_list

    @staticmethod
    def mk_master():
        """ This function will make the master frame that will be used for the differencing

        :return - The master frame is returned and written to the master directory
        """

        file_name = Configuration.FIELD +'_master' + Configuration.FILE_EXTENSION

        if os.path.isfile(Configuration.MASTER_DIRECTORY + file_name) == 0:

            chk_tmp_files = Utils.get_file_list(Configuration.MASTER_TMP_DIRECTORY, Configuration.FILE_EXTENSION)

            # get the image list
            full_image_list, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY,
                                                              Configuration.FIELD,
                                                              'clean')
            # determine the number of loops we need to move through for each image
            full_nfiles = len(full_image_list)

            # get the nights the images were observed
            img_by_night = np.array([line.split('/')[7] for line in full_image_list]).reshape(-1)

            # loop through the image headers and get the sky background
            sky_values = np.zeros(full_nfiles)
            bd_wcs = np.zeros(full_nfiles)
            for idx, file in enumerate(full_image_list):

                # pull in the header file
                h_chk = fits.getheader(file)

                # get the sky background
                sky_values[idx] = h_chk['sky']

                # update the bad index if it exists
                if h_chk['BAD_WCS'] == 'Y':
                        bd_wcs[idx] = 1
                else:
                    bd_wcs[idx] = 0

            # get the statistics on the images
            img_mn, img_mdn, img_std = sigma_clipped_stats(sky_values, sigma=2)

            # now make a list of the images to use for the master frame
            image_list = []
            for idx, file in enumerate(full_image_list):

                # remove nights which have high sky background, bad wcs, or are on hand-picked bad nights
                if (sky_values[idx] <= img_mdn + 2 * img_std) & (bd_wcs[idx] == 0) & (img_by_night[idx] != '2024-10-03'):
                    image_list.append(file)

            nfiles = len(image_list)

            if len(chk_tmp_files) == 0:
                Utils.log("No temporary Master files found. Generating new ones.", "info")
                # determine the number of loops we need to move through for each image
                nbulk = 20

                # get the integer and remainder for the combination
                full_bulk = nfiles // nbulk
                part_bulk = nfiles % nbulk

                if part_bulk > 0:
                    hold_bulk = full_bulk + 1
                else:
                    hold_bulk = full_bulk

                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_Y, Configuration.AXS_X))

                # update the log
                Utils.log("Generating a master frame from multiple files in bulks of " + str(nbulk) +
                          " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                          str(hold_bulk) + " mini-files to median combine.", "info")

                cnt_img = 0
                for kk in range(0, hold_bulk):

                    # loop through the images in sets of nbulk
                    if kk < full_bulk:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(nbulk, Configuration.AXS_Y, Configuration.AXS_X))

                        # generate the max index
                        mx_index = nbulk
                    else:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS_Y, Configuration.AXS_X))

                        # generate the max index
                        mx_index = part_bulk

                    # make the starting index
                    loop_start = kk * nbulk
                    idx_cnt = 0

                    Utils.log("Making mini file " + str(kk) + ".", "info")

                    # now loop through the images
                    for jj in range(loop_start, mx_index + loop_start):
                        # read in the image directly into the block_hold

                        master_tmp, master_tmp_head = fits.getdata(image_list[jj], header=True)

                        if (kk == 0) & (jj == 0):
                            block_hold[idx_cnt] = master_tmp - master_tmp_head['sky']
                            master_header = master_tmp_head
                            del master_tmp
                            Utils.log("Mini file " + str(kk) + " image " + str(jj) +
                                      " aligned. " + str(nfiles - cnt_img) + " remain.",
                                      "info")
                        else:
                            tmp = Preprocessing.align_img(master_tmp, master_tmp_head, master_header)
                            Utils.log("Mini file " + str(kk) + " image " + str(jj) +
                                      " aligned. " + str(nfiles - cnt_img) + " remain.",
                                      "info")
                            block_hold[idx_cnt] = tmp - master_tmp_head['sky']
                            del tmp
                            del master_tmp

                        cnt_img = cnt_img + 1
                        # increase the iteration
                        idx_cnt += 1
                    # median the data into a single file
                    hold_data[kk] = np.median(block_hold, axis=0)
                    del block_hold
                    if kk < 10:
                        fits.writeto(Configuration.MASTER_TMP_DIRECTORY + "0" + str(kk) + "_tmp_master.fits",
                                     hold_data[kk], master_header, overwrite=True)
                    else:
                        fits.writeto(Configuration.MASTER_TMP_DIRECTORY + str(kk) + "_tmp_master.fits",
                                     hold_data[kk], master_header, overwrite=True)
            else:
                Utils.log("Legacy files found. Creating Master frame from these files. "
                          "Delete if you do not want this!", "info")
                hold_bulk = len(chk_tmp_files)
                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_Y, Configuration.AXS_X))

                for kk, tmp_file in enumerate(chk_tmp_files):
                    master_tmp, master_tmp_head = fits.getdata(Configuration.MASTER_TMP_DIRECTORY + tmp_file, header=True)
                    hold_data[kk] = master_tmp
                    master_header = master_tmp_head

            # median the mini-images into one large image
            master = np.median(hold_data, axis=0)

            master_header['MAST_COMB'] = 'median'
            master_header['NUM_MAST'] = nfiles

            # now mask the bad parts of the image #### THIS WILL CHANGE PER FIELD!!!! LIKELY YOU SHOULD REMOVE#####
            master[0:490, :] = 0
            master[10045:-1, :] = 0
            master[:, 0:530] = 0
            master[:, 10465:-1] = 0
            #### END LIKELY REMOVE

            # write the image out to the master directory
            fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                         master, master_header, overwrite=True)
        else:
            master, master_header = fits.getdata(Configuration.MASTER_DIRECTORY + file_name, header=True)

        return master, master_header
