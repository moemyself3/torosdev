import pandas as pd
from config import Configuration
from astropy.io import fits
import numpy as np
from libraries.utils import Utils
biases = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'bias_list.csv', sep=',')
darks = pd.read_csv(Configuration.CALIBRATION_DIRECTORY + 'dark_list.csv', sep=',')

# group based on the date to get the unique dates
bias_dates = biases.groupby('Date').agg({'count'}).reset_index()['Date'].to_list()

f = open('bias_dark_chk.txt', 'w')
f.write('date bias_mean bias_median dark_mean dark_median\n')
for dte in bias_dates:
    bias_frames = biases[biases.Date == dte]
    dark_frames = darks[darks.Date == dte]

    if len(dark_frames) > 0:
        bias_list = bias_frames.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.File, axis=1).to_list()
        mn_bias = np.zeros(len(bias_list))
        mdn_bias = np.zeros(len(bias_list))
        for ii in range(0, len(bias_list)):
            bias_tmp = fits.getdata(bias_list[ii]).astype('float')
            mn_bias[ii] = np.mean(bias_tmp)
            mdn_bias[ii] = np.median(bias_tmp)
        Utils.log('The mean of the bias on ' + dte + ' was ' + str(np.mean(mn_bias)) + '.', "info")
        dark_list = dark_frames.apply(lambda x: Configuration.RAW_DIRECTORY + x.Date + x.Files, axis=1).to_list()
        mn_dark = np.zeros(len(dark_list))
        mdn_dark = np.zeros(len(dark_list))
        for ii in range(0, len(dark_list)):
            dark_tmp = fits.getdata(dark_list[ii]).astype('float')
            mn_dark[ii] = np.mean(dark_tmp)
            mdn_dark[ii] = np.median(dark_tmp)
        Utils.log('The mean of the dark on ' + dte + ' was ' + str(np.mean(mn_dark)) + '.', "info")

        f.write(dte + ' ' + str(int(np.around(np.mean(mn_bias), decimals=0)))
                + ' ' + str(int(np.around(np.median(mdn_bias), decimals=0)))
                + ' ' + str(int(np.around(np.mean(mn_dark), decimals=0)))
                + ' ' + str(int(np.around(np.median(mdn_dark), decimals=0))) + '\n')
f.close()
