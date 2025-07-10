""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats
import numpy as np
import pandas as pd
import warnings
from astropy.time import Time
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
import logging
matplotlib.set_loglevel(level = 'warning')
matplotlib.use("TkAgg")
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.INFO)
import matplotlib.pyplot as plt


class Photometry:

    @staticmethod
    def single_frame_aperture_photometry(star_list, img_name, fin_name):
        """ This function will find the subtraction stars to use for the differencing, they will be the same stars for
        every frame. This will help in detrending later.

        :parameter star_list - The data frame with the list of stars to use for the subtraction
        :parameter img_name - The file to extract flux from
        :parameter fin_name - The name of the output flux file

        :return nothing is returned, but the flux file is written and output
        """

        # get the image for photometry
        img, header = fits.getdata(img_name, header=True)

        # get the various important header information
        time = Time(header['DATE'], format='isot', scale='utc')
        jd = time.jd

        # get the stellar positions from the master frame
        positions = np.transpose((star_list['xcen'], star_list['ycen']))

        aperture = CircularAperture(positions, r=Configuration.APER_SIZE)
        aperture_area = aperture.area
        annulus_aperture = CircularAnnulus(positions,
                                           r_in=Configuration.ANNULI_INNER,
                                           r_out=Configuration.ANNULI_OUTER)

        # get the background stats
        aperstats = ApertureStats(img, annulus_aperture)
        bkg_mean = aperstats.mean
        total_bkg = bkg_mean * aperture_area

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, aperture, method='exact')

        # extract the flux from the table
        # the sky was subtracted during the calibration and differencing steps, the raw photometry should be fine
        star_flux = np.array(phot_table['aperture_sum']) * Configuration.GAIN

        # calculate the expected photometric error
        star_error = np.abs(star_flux)
        bkg_error = np.sqrt((header['sky'] * aperture_area) ** 2 + total_bkg ** 2) * Configuration.GAIN

        # combine sky and signal error in quadrature
        star_flux_err = np.sqrt(star_error + bkg_error)
        star_flux_snr = np.abs(star_flux / star_flux_err)

        # combine the fluxes
        flux = star_flux.astype(float) + star_list['master_flux'].to_numpy().astype(float)
        flux_er = np.sqrt(star_flux_err.astype(float) ** 2 + star_list['master_flux_er'].to_numpy().astype(float) ** 2)

        # convert to magnitude
        mag = 25. - 2.5 * np.log10(flux)
        mag_er = (np.log(10.) / 2.5) * (flux_er / flux)

        # get the zeropoint
        m_mag = star_list['master_mag'].to_numpy()

        dmag = mag[~np.isnan(mag)] - m_mag[~np.isnan(mag)]

        f_mags = np.arange(6, 16) + 0.5
        n_mags = np.zeros(len(f_mags))
        for mag_idx, m_lw in enumerate(f_mags):
            dmag_bin = dmag[np.argwhere((mag[~np.isnan(mag)] > m_lw) & (mag[~np.isnan(mag)] < m_lw + 1))]
            dmag_mn, dmag_md, dmag_sg = sigma_clipped_stats(np.array(dmag_bin, dtype=float), sigma=2.5)
            n_mags[mag_idx] = dmag_md

        # replace nans with -9.999999 and calculate the offset
        mag = np.where(np.isnan(mag), -9.999999, mag)
        off = np.interp(mag, f_mags, n_mags)

        # now correct the magnitudes for exposure time
        mag = mag + 2.5 * np.log10(Configuration.EXP_TIME)

        # generate the final flux file
        flux_file = star_list.copy().reset_index(drop=True)
        flux_file['flux'] = flux
        flux_file['flux_er'] = flux_er
        flux_file['mag'] = mag
        flux_file['mag_er'] = mag_er
        flux_file['sky'] = header['SKY']
        flux_file['jd'] = jd
        flux_file['zpt'] = off
        flux_file['exp_time'] = Configuration.EXP_TIME

        flux_file.to_csv(fin_name, header=True, index=False)
        return

    @staticmethod
    def combine_flux_files():
        """ This function combines all of the flux files in a given directory into a single data frame.

        :return - Nothing is being returned, but the raw files are output to disk
        """

        # pull in the star list for the photometry
        star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt',
                                delimiter=' ',
                                header=0)

        # get the flux files to read in
        files, dates = Utils.get_all_files_per_field(Configuration.DIFFERENCED_DIRECTORY,
                                                     Configuration.FIELD,
                                                     'diff',
                                                     '.flux')
        nfiles = len(files)

        num_rrows = len(star_list)

        # make the holders for the light curves
        jd = np.zeros(nfiles)
        mag = np.zeros((num_rrows, nfiles))
        er = np.zeros((num_rrows, nfiles))
        trd = np.zeros((num_rrows, nfiles))
        zpt = np.zeros((num_rrows, nfiles))

        for idy, file in enumerate(files):

            # read in the data frame with the flux information
            img_flux = pd.read_csv(file, header=0)

            if idy == 0:
                src_id = img_flux['source_id'].to_numpy()

            # set the data to the numpy array
            jd[idy] = img_flux.loc[0, 'jd']
            mag[:, idy] = img_flux['mag'].to_numpy()
            er[:, idy] = img_flux['mag_er'].to_numpy()
            zpt[:, idy] = img_flux['zpt'].to_numpy()

            if (idy % 100 == 0) & (idy > 0):
                Utils.log("100 flux files read. " + str(nfiles - idy - 1) + ' files remain.', "info")

        for idy, row in star_list.loc[0:30000].iterrows():

            # get the distance to all stars
            dd = np.sqrt((row.xcen - star_list.xcen.to_numpy()) ** 2 +
                         (row.ycen - star_list.ycen.to_numpy()) ** 2)

            # get the difference in magnitude
            dmag = np.abs(row.master_mag - star_list.master_mag.to_numpy())

            # only get nearby stars of similar magnitude
            vv = np.argwhere((dd < 500) & (dd > 0) & (dmag > 0) & (dmag < .5)).reshape(-1)

            if len(vv) > 0:
                # make the trend holder and standard deviation holder
                holder = np.zeros((len(vv), len(jd)))

                # loop through all OK stars removing outliers
                for idz in range(len(vv)):
                    _, mdn, _ = sigma_clipped_stats(mag[vv[idz], :], sigma=2.5)
                    holder[idz, :] = mag[vv[idz], :] - mdn

                for idz in range(len(jd)):
                    _, trd[idy, idz], _ = sigma_clipped_stats(holder[:, idz], sigma=2)

            if (idy % 1000 == 0) & (idy > 0):
                Utils.log("1000 stars had their trend found. " +
                          str(num_rrows - idy - 1) + " stars remain.", "info")

        # write out the light curve data
        Photometry.write_light_curves(num_rrows, jd, mag, er, trd, zpt, src_id)

        return

    @staticmethod
    def write_light_curves(nstars, jd, mag, er, trd, zpt, src_id):
        """ This function will write the flux columns to light curves for each src_id

        :parameters - nstars - the number of stars to make light curves
        :parameters - jd - the numpy array of julian dates (one per file)
        :parameters - mag - the magnitudes for each star at each jd
        :parameters - er - the photometric error for each star at each time
        :parameters - trd - the trend from nearby similar magnitude stars
        :parameter - zpt - the zeropoint from all photometry
        :paramters - src_id - nstars long array of source ids

        :return - Nothing is returned, but the light curve files are written
        """

        # initialize the light curve data frame
        lc = pd.DataFrame(columns=['jd', 'mag', 'er', 'cln', 'zpt'])

        Utils.log("Starting light curve writing...", "info")

        for idx in range(0, 30000):
            star_id = str(src_id[idx])

            # add the time, magnitude and error to the data frame
            lc['jd'] = np.around(jd, decimals=6)
            lc['mag'] = np.around(mag[idx, :], decimals=6)
            lc['er'] = np.around(er[idx, :], decimals=6)
            lc['trd'] = np.around(trd[idx, :], decimals=6)
            lc['zpt'] = np.around(zpt[idx, :], decimals=6)

            lc = lc.sort_values(by = 'jd')
            # write the new file
            lc[['jd', 'mag', 'er', 'trd', 'zpt']].to_csv(Configuration.LIGHTCURVE_DIRECTORY +
                                                         Configuration.FIELD + "/" +
                                                         Configuration.FIELD +"_" + star_id + ".lc",
                                                         sep=" ", index=False, na_rep='9.999999')

        return
