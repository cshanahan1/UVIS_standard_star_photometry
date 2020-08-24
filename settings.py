"""
Configuration file for standard_star_phot_pipeline. Set all paths and parameters
here before running standard_star_phot_pileline.py.
"""

import os

# Path to directory where data will be downloaded, and sorted into 
DATA_DIR = '/grp/hst/wfc3p/cshanahan/phot_group_work/data/standard_star_data/staring_mode_data'
DATA_DIR = os.path.join(DATA_DIR, '')

# output directory for log files, photometry catalots
OUTPUT_DIR = '/grp/hst/wfc3p/cshanahan/phot_group_work/analysis_pipelines/standard_star_photometry_pipeline/output/'
OUTPUT_DIR = os.path.join(OUTPUT_DIR, '')

# Desired file type for download/photometry (i.e flc, flt, drz...)
FILE_TYPE = 'flc'

# Set list of proposal IDs. Used for MAST query when downloading data.
PROPOSAL_IDS = [11426, 11450, 11557, 11903, 11907, 12090, 
				12333, 12334, 12698, 12699, 12707, 13088, 
				13089, 13096, 13574, 13575, 13584, 13711, 
				14018, 14021, 14382, 14384, 14815, 14883, 
				14992, 15113, 15398, 15399, 15582, 15583]

# Set dictionary of all target names, and possible variations.
TARGNAME_MAPPINGS = {'G191B2B' : ['G191B2B'],
					 'GD153' : ['GD153', 'GD-153'],
					 'GRW70' : ['GRW70', 'GRW+70D5824', 'GRW+70D'],
					 'GD71' : ['GD71', 'GD-71'],
					 'P330E' : ['P330E', 'GSC-02581-02323']}
ALL_TARGETS = [i for x in [v for k, v in TARGNAME_MAPPINGS.items()] for i in x]

# Boolean, if LACosmic should be run and CRs masked in data before photometry
RUN_CR_REJECTION = False

# Set list of desired circular aperture radii for photometry
CIRCULAR_APERTURE_RADII = [3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 
						   9.5, 10, 12, 15, 17, 20, 25, 30, 40, 50, 60, 70, 80, 
						   90, 100]

# Set (inner sky annulus radius, outer sky annulus radius) for sky region
SKY_ANNULUS = (156., 165.)

# Choice of 'mean', 'median', 'mode' for computing sky level. 
SKY_SUB_METHOD = 'mode'

# Boolean, if sky region should be sigma clipped
SIGMA_CLIP_SKY = True

# Directory containing pixel area maps
PAM_PATH = '/grp/hst/wfc3p/cshanahan/phot_group_work/pixel_area_maps/'