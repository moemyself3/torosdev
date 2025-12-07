#This code is adated from Reyes, A. (2024). On the Development of a Software Pipeline for Automatic Detection of Transients With the Transient Optical Robotic Observatory of the South [Master's thesis, The University of Texas Rio Grande Valley]. ScholarWorks @ UTRGV. https://scholarworks.utrgv.edu/etd/1607

from config import Configuration
from libraries.utils import Utils
import os
import copy
from math import ceil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

from astropy.io import fits
from astropy.table import Table, setdiff
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from photutils.detection import find_peaks, DAOStarFinder
from photutils.psf import extract_stars, EPSFBuilder
from photutils.psf.matching import resize_psf
from photutils.aperture import ApertureStats, CircularAperture
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree

def patch(position, image, fake):

    position = (position[1], position[0]) # first element is y below
    coords = [(int(position[0]-fake.shape[0]/2),int(position[0]+fake.shape[0]/2)),
              (int(position[1]-fake.shape[1]/2),int(position[1]+fake.shape[1]/2))]

    patch = np.add(image[coords[0][0]:coords[0][1],coords[1][0]:coords[1][1]],fake)
    image[coords[0][0]:coords[0][1],coords[1][0]:coords[1][1]] = patch

    return image

def write_fits(data, header, filename, overwrite=False):
    hdul = fits.HDUList()
    hdul.append(fits.ImageHDU(data=data, header=header))
    hdul.writeto(filename, overwrite=overwrite)

class Injector:
    """Simulated source injection class.
        Select candidate stars with find_peaks() for PSF modelling
        using photutils.psf.EPSFBuilder().
        Re-scales the PSF model and injects at specified location and magnitude.

    Arguments:
    image : `str`
        path to the fits image subject to injection
    threshold : `int`, default=500
        photutils find_peaks() detection threshold
    cutout_size : `int`, default=30
        size of candidate star square cutout sides

    Attributes:
    -----------
    epsf_resized : `numpy.ndarray`
        PSF model resized to match cutout_size
    """

    def __init__(self, image, threshold=100, 
                 cutout_size = None,
                 epsf_resized = None
                 ) -> None:

        self.image = image

        # load image
        try:
            self.hdul = fits.open(self.image, mode='readonly')
            for i in range(10):
                if type(self.hdul[i].data) != type(None):
                    self.image = self.hdul[i]
                    break
        except FileNotFoundError as e:
            Utils.log(f"{e}", "error")

        self.threshold = threshold

        if cutout_size is None:
            cutout_size = ceil(Configuration.APER_SIZE*2.5)

        if cutout_size % 2 == 0:
            # if even
            cutout_size += 1 # make odd

        self.cutout_size = cutout_size
        self.epsf = None
        self.epsf_resized = epsf_resized
        self.fitted_stars = None
        self.lonely_peaks_table = None
        self.image_median = None

    def make_candidate_starlist(self):
        """Generates a candidate list of stars to use to make epsf.

        Arguments:
        ----------
        None

        Returns:
        --------
        lonely_peaks_table

        """

        # Extract sources x y values
        peaks_DAOStarFinder = DAOStarFinder(
                                fwhm=Configuration.FWHM,
                                threshold=self.threshold,
                                exclude_border=True,
                                peakmax=Configuration.PEAKMAX)
        weaks_DAOStarFinder = DAOStarFinder(
                                fwhm=Configuration.FWHM, 
                                threshold=self.threshold*0.4,
                                exclude_border=True,
                                peakmax=Configuration.PEAKMAX)
        Utils.log("start sigma_clipped_stats", "debug")
        mean, median, std = sigma_clipped_stats(self.image.data,
                                                cenfunc='median',
                                                stdfunc='std')
        self.image_median = median
        median_subtracted_img = self.image.data - median
        Utils.log("Start peaks weaks DAOStarFinder", "debug")
        peaks_table = peaks_DAOStarFinder(median_subtracted_img)
        weaks_table = weaks_DAOStarFinder(median_subtracted_img)
        Utils.log("End peaks weaks DAOStarFinder", "debug")

        # keep only "real" stars
        peaks_table = peaks_table[(peaks_table['flux'] > 0) ]#& (peaks_table['sharpness'] < 0.7)]
        weaks_table = weaks_table[(weaks_table['flux'] > 0) ]#& (weaks_table['sharpness'] < 0.7)]
        peaks_table.rename_column('xcentroid', 'x_peak')
        peaks_table.rename_column('ycentroid', 'y_peak')
        weaks_table.rename_column('xcentroid', 'x_peak')
        weaks_table.rename_column('ycentroid', 'y_peak')
        peaks_table['x_peak'] = peaks_table['x_peak'].astype(int)
        peaks_table['y_peak'] = peaks_table['y_peak'].astype(int)
        weaks_table['x_peak'] = weaks_table['x_peak'].astype(int)
        weaks_table['y_peak'] = weaks_table['y_peak'].astype(int)

        # remove peaks_table sources from weaks_table
        weaks_table = setdiff(weaks_table, peaks_table, keys=['x_peak','y_peak'])

        # Remove neighboring sources
        lonely_mask = []
        # get distances from each peaks_table source to weaks_table source
        peaks_xy_array = peaks_table['x_peak', 'y_peak'].to_pandas().to_numpy().astype(int)
        weaks_xy_array = weaks_table['x_peak', 'y_peak'].to_pandas().to_numpy().astype(int)

        chunk_size = 10000
        nchunks = len(peaks_xy_array) // chunk_size
        if nchunks == 0:
            nchunks = 1

        peaks_xy_chunks =np.array_split(peaks_xy_array, nchunks)
        for idx, peaks_xy_chunk in enumerate(peaks_xy_chunks):
            if idx == 0:
                distances = cdist(peaks_xy_chunk, weaks_xy_array)
            else:
                Utils.log(f'chunking peaks distance calculation: {idx}', 'info')
                np.append(distances, cdist(peaks_xy_chunk, weaks_xy_array))

        peaks_min_distances = np.min(distances, axis=1)
        peaks_table['min_distances'] = peaks_min_distances
        lonely_mask = peaks_min_distances > self.cutout_size

        lonely_peaks_table = peaks_table[lonely_mask]

        del peaks_table, weaks_table, peaks_xy_array, weaks_xy_array, distances, peaks_min_distances

        # remove lonely peaks less than the cutout_size to avoid contamination
        lonely_xy_array = lonely_peaks_table['x_peak', 'y_peak'].to_pandas().to_numpy().astype(int)
        tree = KDTree(lonely_xy_array)
        distances, indices = tree.query(lonely_xy_array, k=2)
        min_distances_kdtree = distances[:, 1]
        lonely_peaks_table =  lonely_peaks_table[ (min_distances_kdtree > self.cutout_size) ]
        # add image_median into metadata of table
        lonely_peaks_table.meta['image_median'] = self.image_median
        self.lonely_peaks_table = lonely_peaks_table
        return lonely_peaks_table

    def make_epsf(self, psf_stats=True, verb=True):
        """Generates effective PSF from image using photutils.psf.EPSFBuilder()

        Arguments:
        ----------
        psf_stats : `bool`, default=True
            computes the magnitude of the psf model using photutils.psf.ApertureStats() for later injection
        """
        Utils.log("Making EPSF.", "info")
        if self.lonely_peaks_table is None:
            lonely_peaks_table = self.make_candidate_starlist()
        else:
            lonely_peaks_table = self.lonely_peaks_table

        # Vignet mask
        hcutout_size = int(self.cutout_size/2)
        shape_x = self.image.data.shape[1]
        shape_y = self.image.data.shape[0]
        x = lonely_peaks_table['x_peak']
        y = lonely_peaks_table['y_peak']
        mask = ((x > hcutout_size + 0.1*shape_x) & (x < (0.9*shape_x -1 - hcutout_size)) &
                (y > hcutout_size + 0.1*shape_y) & (y < (0.9*shape_y -1 - hcutout_size)))

        # Good star positions table
        star_table = Table()
        star_table['x'] = x[mask]
        star_table['y'] = y[mask]

        # Bad pixel mask using off-center pixel
        bp_mask = []
        for i in range(len(star_table['x'])):
            bp_mask.append(True)
            from_center = int(self.cutout_size/6) # pixel 1/3 from bottom left corner
            thresh_val = self.threshold * 0.4
            if self.image.data[star_table['y'][i]-from_center][star_table['x'][i]-from_center] < thresh_val:
                bp_mask[i] = False

        star_table = star_table[bp_mask]

        # Extract stars
        median_subtracted_img = self.image.data - self.image_median
        nddata = NDData(median_subtracted_img)
        stars = extract_stars(nddata, star_table, size=self.cutout_size)

        # Initiate EPSFBuilder
        Utils.log("Starting EPSFBuilder.", "info")
        epsf_builder = EPSFBuilder(oversampling=4, maxiters=5, progress_bar=False,
                                                smoothing_kernel='quadratic')
        self.epsf, self.fitted_stars = epsf_builder(stars)
        self.stars = stars

        self.epsf_resized = resize_psf(self.epsf.data, 1, len(self.epsf.data)/self.cutout_size)
        self.epsf_resized[self.epsf_resized < 0] = 0

        if psf_stats == True:
            Utils.log("Calculating PSF stats", "info")
            psf_side = self.epsf_resized.shape[0]
            psf_aperture = CircularAperture((int(psf_side/2),int(psf_side/2)), r=psf_side*.5)
            psf_stats = ApertureStats(self.epsf_resized.data, psf_aperture)
            psf_fwhm = psf_stats.fwhm.to_value()
            # Improved radius using fwhm:
            self.psf_aperture = CircularAperture((int(psf_side/2),int(psf_side/2)), r=psf_fwhm*1.5)
            self.psf_stats = ApertureStats(self.epsf_resized.data, self.psf_aperture)
            #self.psf_magnitude = -2.5*np.log10(self.psf_stats.sum/self.image.header['EXPTIME'])
            self.psf_magnitude = -2.5*np.log10(self.psf_stats.sum)

    def load_epsf(self, epsf_resized):
            self.epsf_resized = epsf_resized
            psf_side = self.epsf_resized.shape[0]
            psf_aperture = CircularAperture((int(psf_side/2),int(psf_side/2)), r=psf_side*.5)
            psf_stats = ApertureStats(self.epsf_resized.data, psf_aperture)
            psf_fwhm = psf_stats.fwhm.to_value()
            self.psf_aperture = CircularAperture((int(psf_side/2),int(psf_side/2)), r=psf_fwhm*1.5)
            self.psf_stats = ApertureStats(self.epsf_resized.data, self.psf_aperture)
            self.psf_magnitude = -2.5*np.log10(self.psf_stats.sum)

    def inspect_stars(self, title=''):
        """Plot stars used for modeling EPSF

        Arguments:
        ----------
        title : `str`, default=''
            plot title
        """

        if self.fitted_stars == None:
            Utils.log('No sources found. Run make psf first.','info')

        else:
            fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(35,35), squeeze=True)
            ax = ax.ravel()
            ncols, nrows = (5,5)

            for i in range(nrows*ncols):
                j = int(len(self.fitted_stars)*i/len(range(nrows*ncols)))
                norm = simple_norm(self.fitted_stars[j].data, 'log', percent=99.)
                ax[i].imshow(self.fitted_stars[j].data, norm=norm, origin='lower', cmap='viridis')

            if title == '':
                title = 'Candidate stars for PSF modelling (25 of %i)' % int(len(self.fitted_stars))
            plt.suptitle(title)
            plt.show()

    def plot_psf(self, title=''):
        """Plots 3D projection of PSF model.
        """

        psf_asarray = np.asarray(self.epsf_resized.data)
        x = np.arange(0,self.cutout_size)
        y = np.arange(0,self.cutout_size)
        x,y = np.meshgrid(x,y)

        fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
        surf = ax.plot_surface(x, y, psf_asarray, cmap=cm.coolwarm, antialiased=False)

        ax.set_zlim(0,1.02*psf_asarray.max())
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter('{x:.01f}')
        ax.set_title(title)
        ax.set_xlabel('x axis')
        ax.set_ylabel('y axis')
        ax.set_zlabel('counts')

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()

        return

    def inject(self, position, magnitude, image):#=None):
        """Inject simulated source into location in image

        Arguments:
        ----------
        position : `tuple`
            (x, y) `int` pixel position of simulated source
        magnitude: `int`, `float`
            target magnitude of simulated source
        image : `numpy.ndarray`, default=None
            image where simulated source will be injected
            if None, copy of original image will be used
        return_data: `bool`

        Returns:
        --------
        image : `np.array`, optional
            image with injected simulated source
        """

        if type(image) == type(None):
            image = copy.deepcopy(self.image.data)
            Utils.log('Using input fits image data', 'info')
        # else:
        #     image = image

        #target_counts = self.image.header['EXPTIME']*10**(-magnitude/2.5)
        target_counts = 10**(-magnitude/2.5)
        scalar = target_counts / self.psf_stats.sum

        sim_source = self.epsf_resized.data * scalar
        return patch(position, image, sim_source)
