import pandas as pd
import matplotlib
import logging
matplotlib.set_loglevel(level = 'warning')
matplotlib.use("TkAgg")
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.INFO)
import matplotlib.pyplot as plt
from config import Configuration
import numpy as np
from libraries.utils import Utils
from astropy.stats import sigma_clipped_stats as scs

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
trd = np.zeros(nfiles)
zpt = np.zeros((num_rrows,nfiles))

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

idy = 25706

# get the distance to all stars
dd = np.sqrt((star_list.xcen[idy] - star_list.xcen.to_numpy()) ** 2 +
             (star_list.ycen[idy] - star_list.ycen.to_numpy()) ** 2)

# get the difference in magnitude
dmag = np.abs(star_list.master_mag[idy] - star_list.master_mag.to_numpy())

# only get nearby stars of similar magnitude
vv = np.argwhere((dd < 500) & (dd > 0) & (dmag > 0) & (dmag < .5)).reshape(-1)

# make the trend holder and standard deviation holder
holder = np.zeros((len(vv), len(jd)))

# loop through all OK stars removing outliers
for idz in range(len(vv)):
    _, mdn, _ = scs(mag[vv[idz], :], sigma=2.5)
    holder[idz, :] = mag[vv[idz], :] - mdn

for idz in range(len(jd)):
    _, trd[idz], _ = scs(holder[:, idz], sigma=2)


raw = mag[idy,:]
ers = er[idy,:]
img_list = pd.read_csv(Configuration.CLEAN_DIRECTORY + "image_list.csv", sep=',')

bd = img_list['Bad'].to_numpy()
plt.scatter(jd[(bd == 0) & (mag > 0)] - 2460584, raw[(bd == 0) & (mag > 0)] - trd[(bd == 0) & (mag > 0)], marker='.', c='k')
plt.scatter(jd[bd == 1] - 2460584, raw[bd == 1] - trd[bd == 1], marker='.', c='r')
plt.xlabel('JD - 2460584 [d]')
plt.ylabel('T')

plt.gca().invert_yaxis()
plt.show()


star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.FIELD + '_star_list.txt',
                        delimiter=' ',
                        header=0)

lc = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.FIELD + "/" + Configuration.FIELD +"_4689579240377005184.lc", sep=' ')
lc['ph'] = lc.apply(lambda x: (x.jd - np.min(lc.jd.to_numpy())) / 0.37143 % 1, axis=1)
lc['bd'] = img_list['Bad']
lc['bd'] = np.where(lc.mag < 0, 1, lc.bd)

plt.subplot(2, 1, 1)
plt.errorbar(lc[lc.bd == 0].ph, lc[lc.bd == 0].mag - lc[lc.bd == 0].trd, yerr=lc[lc.bd == 0].er, c='k', fmt='None')
plt.scatter(lc[lc.bd == 0].ph, lc[lc.bd == 0].mag - lc[lc.bd == 0].trd, marker='.', c='k')
plt.scatter(lc[lc.bd == 1].ph, lc[lc.bd == 1].mag - lc[lc.bd == 1].trd, marker='.', c='r')
plt.xlabel('Phase')
plt.ylabel('T')
plt.title('CO Tuc P = 0.37143 d')
plt.gca().invert_yaxis()

plt.subplot(2, 1, 2)
plt.errorbar(lc[lc.bd == 0].jd - 2460584, lc[lc.bd == 0].mag - lc[lc.bd == 0].trd, yerr=lc[lc.bd == 0].er, c='k', fmt='None')
plt.scatter(lc[lc.bd == 0].jd - 2460584, lc[lc.bd == 0].mag - lc[lc.bd == 0].trd, marker='.', c='k')
plt.scatter(lc[lc.bd == 1].jd - 2460584, lc[lc.bd == 1].mag - lc[lc.bd == 1].trd, marker='.', c='r')
plt.xlabel('JD - 2460584 [d]')
plt.ylabel('T')
plt.gca().invert_yaxis()
plt.show()

lcs = Utils.get_file_list(Configuration.LIGHTCURVE_DIRECTORY + Configuration.FIELD + "/", '.lc')

for idx, ll in enumerate(lcs):

    if idx % 100 == 0:
        print(idx)
        lc = pd.read_csv(
            Configuration.LIGHTCURVE_DIRECTORY + Configuration.FIELD + "/" + ll,
            sep=' ')
        lc['bd'] = img_list['Bad']
        lc['bd'] = np.where(lc.mag < 0, 1, lc.bd)

        _, mdn, sts = scs(lc[lc.bd == 0].mag - lc[lc.bd == 0].zpt, sigma=2.5)
        plt.errorbar(lc[lc.bd == 0].jd - 2460584, lc[lc.bd == 0].mag - lc[lc.bd == 0].zpt, yerr=lc[lc.bd == 0].er, c='k', fmt='None')
        plt.scatter(lc[lc.bd == 0].jd - 2460584, lc[lc.bd == 0].mag - lc[lc.bd == 0].zpt, marker='.', c='k')
        plt.scatter(lc[lc.bd == 1].jd - 2460584, lc[lc.bd == 1].mag - lc[lc.bd == 1].zpt, marker='.', c='r')
        plt.xlabel('JD - 2460584 [d]')
        plt.ylabel('T')
        plt.ylim([mdn + 2.5 * sts, mdn - 2.5 * sts])
        plt.show()