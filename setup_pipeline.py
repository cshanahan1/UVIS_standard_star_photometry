import os

from settings import *

""" Sets up output and data directories/subdirectories for photometry pipeline. 
	Creates DATA_DIR, set in settings.py, and other subdirectories if they don't
	exist. """
	
def setup_pipeline():

	dirs_to_make = [DATA_DIR, OUTPUT_DIR, OUTPUT_DIR+'photcats']
	for d in dirs_to_make:
		if not os.path.isdir(d):
			print('Making {}.'.format(d))
			os.makedirs(d)

if __name__ == '__main__':
	setup_pipeline()