import pandas as pd
from astropy.io import fits
import numpy as np
from libraries.utils import Utils
from config import Configuration
from astropy.time import Time
from astropy.stats import sigma_clip, mad_std
_, dates = Utils.get_all_files_per_field(Configuration.CLEAN_DIRECTORY,
                                         Configuration.FIELD,
                                         Configuration.FILE_EXTENSION)

for dte in dates:
    cln_dir = Configuration.CLEAN_DIRECTORY + dte + '/' + Configuration.FIELD +'/'
    wcs_df = pd.DataFrame(columns=['JD', 'CRPIX1', 'CRPIX2', 'PC1_1', 'PC1_2',
                                   'PC2_1', 'PC2_2', 'CRVAL1', 'CRVAL2', 'LATPOLE'])
    all_files = Utils.get_file_list(cln_dir, Configuration.FILE_EXTENSION)
    nfiles = len(all_files)

    if len(all_files) > 0:
        all_files.sort()

        for idx, file in enumerate(all_files):
            chk_header = fits.getheader(cln_dir + file)

            time = Time(chk_header['DATE'], format='isot', scale='utc')
            jd = time.jd
            try:
                if chk_header['BAD_WCS'] == 'Y':
                    continue
            except:
                wcs_df.loc[idx] = np.array([jd, chk_header['CRPIX1'], chk_header['CRPIX2'], chk_header['PC1_1'],
                                            chk_header['PC1_2'], chk_header['PC2_1'], chk_header['PC2_2'],
                                            chk_header['CRVAL1'], chk_header['CRVAL2'], chk_header['LATPOLE']])
        fd = wcs_df.index.to_numpy()
        fd1 = sigma_clip(wcs_df.CRPIX1.to_numpy(), sigma=3, stdfunc=mad_std).mask
        fd2 = sigma_clip(wcs_df.CRPIX2.to_numpy(), sigma=3, stdfunc=mad_std).mask
        fd3 = sigma_clip(wcs_df.CRVAL1.to_numpy(), sigma=3, stdfunc=mad_std).mask
        fd4 = sigma_clip(wcs_df.CRVAL2.to_numpy(), sigma=3, stdfunc=mad_std).mask

        bd_idx = fd[np.argwhere((fd1 == True) | (fd2 == True) | (fd3 == True) | (fd4 == True))]

        if len(bd_idx) > 0:

            for bd in bd_idx:
                wcs_df = wcs_df.drop(index=bd[0])

            full_idx = np.arange(nfiles)

            for idy in full_idx:
                if idy in wcs_df.index:
                    continue
                else:
                    try:

                        chk_img, chk_header = fits.getdata(cln_dir + all_files[idy], header=True)
                        if chk_header['BAD_WCS'] == 'Y':
                            og_header = chk_header.copy()
                            time = Time(chk_header['DATE'], format='isot', scale='utc')
                            jd = time.jd

                            header_cards = ['CRPIX1', 'CRPIX2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CRVAL1', 'CRVAL2', 'LATPOLE']

                            for hc in header_cards:
                                vv = np.polyfit(wcs_df['JD'], wcs_df[hc], 1)
                                pp = np.poly1d(vv)
                                chk_header[hc] = pp(jd)

                            fits.writeto(Configuration.DATA_DIRECTORY + 'align_test/' + all_files[idy], chk_img, chk_header, overwrite=True)
                    except:
                        continue