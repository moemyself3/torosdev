# Forked toros pipeline

This directory is a development fork of [ryanoelkers/toros](https://github.com/ryanoelkers/toros).

## Major Requirements
- python 3.11
- astrometry net
- source extractor

For a more complete list of requirements, check `toros_env_linux.yml`.


## Quickstart
- Modify config.py to update the location of your data.
    - WORKING\_DIRECTORY
    - DATA\_DIRECTORY
    - RAW\_DIRECTORY
- run `python main.py`

## What's different

### Directory Structure

[ryanoelkers/toros](https://github.com/ryanoelkers/toros) version uses the raw directory directly for raw data, bias, darks, and flats:

 `Configuration.RAW_DIRECTORY + x.Date + x.File`

This version uses the data directory and subfolders for raws, bias, darks, and flats: 

` Configuration.DATA_DIRECTORY + "bias/"  + x.Date + x.File`
` Configuration.DATA_DIRECTORY + "darks/" + x.Date + x.Files`
` Configuration.DATA_DIRECTORY + "flats/" + x.Date + x.File`

#### 

```
RAW_DIRECTORY
├── DATE_1
│   ├── FIELD_1
│   │   ├── FIELD_1_001.fits
│   │   ├── FIELD_1_002.fits
│   │   ├── FIELD_1_003.fits
│   ├── FIELD_2
│   │   ├── FIELD_2_005.fits
│   ├── FIELD_3
│   │   ├── FIELD_3_006.fits
├── DATE_2
│   ├── FIELD_2
│   │   ├── FIELD_2_007.fits
│   ├── FIELD_4
│   │   ├── FIELD_4_008.fits
│   ├── FIELD_5
│   │   ├── FIELD_5_009.fits


DATA_DIRECTORY
├── calibration
│   ├── bias_list.csv
│   ├── dark_list.csv
│   ├── flat_list.csv
```

### Other significant changes
- There was a bug related to the location of the flux files, though that could have been a design change made on the original repo, and I just needed to update my codebase to match.
- modified libraries.utils method `get_all_files_per_field()` to account for the different directory structure. The change in directory structure did affect libraries.preprocessing in a few spots, for example, `mk_flat()` , `mk_combined_bias_and_dark()` , and `mk_nme()` , because the logic used in making the paths was modified.
- opted for astrometry\_net combined with bertin’s source extractor for WCS, in place of twirl.
- modified scripts.difference.BigDiff method `diff_img()` to account for the WCS method used with astrometry\_net (CD\_ij + distortions) vs twirl (PC\_ij). It should work with both.
- modified scripts.difference.BigDiff method `prep_ois()` to work with the linkers called on by gcc. I had issues with how the libraries were called on Mac vs Linux. It should work on both now.
- save to file during master dark, bias, flat field generation instead of keeping chunks in memory.
- modified libraries.preprocessing method `sky_subtract()` by commenting out specific masks that were used for a crowded field specific to a frame.
- added optional parallelization to cleaning and differencing
- added some config variables and directories

## Known Issues
- There are some odd behaviours when using `twirl` to platesolve. Using `astrometry_net` with `source extractor` had better results.
- Upon first install, the working directory needs a `toros.log` file under the `logs` directory. It is not automatically made.
- Inconsistent headers in calibration frame lists. File (in bias\_list.csv and flat\_list.csv) vs Files (in dark\_list.csv)!
