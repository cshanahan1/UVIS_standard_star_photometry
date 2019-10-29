from astropy.io import fits 
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, vstack
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import copy
import glob
import numpy as np
import os 

from settings import *
from WFC3_phot_tools.data_tools import get_wfc3_data_astroquery
from WFC3_phot_tools.data_tools import sort_data
from WFC3_phot_tools.photometry import aperture_phot, utils

def _find_existing_files():
	"""Looks for already downlaoded files in expected subdirectory structure.
	   Returns list of rootnames of files that have already been 
	   downloaded."""

	existing_dat = glob.glob(DATA_DIR+ \
								'*/*/*/*{}.fits'.format(FILE_TYPE))
	# also unsorted data
	existing_dat += glob.glob(DATA_DIR+'*{}.fits'.format(FILE_TYPE))
	# data in 'exclude' directory, don't download 
	existing_dat += glob.glob(DATA_DIR+\
							  'exclude/*{}.fits'.format(FILE_TYPE))
	existing_roots = [os.path.basename(x)[0:9] for x in existing_dat]
	return existing_roots

def _remove_unwanted_files():
	"""Removes all files of in `DATA_DIR` that are 1) spatial scans 
	  2) non-subarray or 3) grism observations, of type `FILE_TYPE`."""

	print('Checking for grism, scan, or non-subarray data.')

	new_data_files = glob.glob(DATA_DIR+'*{}.fits'.format(FILE_TYPE))
	n = 0
	for f in new_data_files:
		hdr = fits.open(f)[0].header
		if hdr['scan_typ'] != 'N':
			print('Removing {}, spatial scan.'.format(f))
			n += 1
			os.remove(f)
			continue
		if 'G' in hdr['filter']:
			print('Removing {}, grism.'.format(f))
			n += 1
			os.remove(f)
			continue
		subarray = hdr['aperture']
		if (subarray == 'UVIS1') or (subarray == 'UVIS2') or (subarray == 'UVIS'):
			print('Removing {}, wrong subarray.'.format(f))
			n += 1
			os.remove(f)
	print('Removed {} files.'.format(n))

def _add_hdr_values_to_phottab(phot_tab, hdr, keywords=['rootname', 'proposid', 
													   'date-obs', 'expstart', 
													   'exptime', 'ccdamp', 
													   'aperture']):
	for keyword in keywords:
		phot_tab[keyword] = hdr[keyword]
	return phot_tab

def _clip_and_log_flux_outliers(phot_tab, threshold, targ, filt, col='countrate_10'):
	"""Clips flux values from phot_tab that are more than `threshold` percent
	higher or lower than the median flux value. 

	Uses the flux in a 10 pixel aperture radius for this calculation by 
	default - for a different flux column, set 'col'. Removes outliers from 
	table and logs them to file in DATA_DIR."""

	median_flux = np.median(phot_tab[col])
	remove_rows = []
	
	with open(OUTPUT_DIR + 'flux_outliers.dat', 'a') as fl:
		for i, row in enumerate(phot_tab):
			percent_dif = np.abs(((row[col]-median_flux)/median_flux)*100)
			if percent_dif >= threshold:
				print('Outlier in {}, {} percent above median.'.format(row['rootname'], percent_dif))
				remove_rows.append(i)
				fl.write('{}, {}, {}, {}, {}, {}\n'.format(row['rootname'], 
					     targ, filt, row['xcenter'], row['ycenter'], 
					     row[col], percent_dif))

	if len(remove_rows) > 0:
		phot_tab.remove_rows(remove_rows)
		print('Removed {} outlier values from photometry table.'.format(len(remove_rows)))

def main_standard_star_phot_pipeline(get_new_data=False, redownload_data=False,
									 sort_new_data=False, cr_reject=False,
									 run_ap_phot=True,
									 pamcorr=True,
									 show_source_detection_plot=False,
									 clip_flux_outliers=True):
	###########################################################################		
	############################ download new data ############################
	###########################################################################	
	if get_new_data:
		query_products = get_wfc3_data_astroquery.query_by_propid_targ_filter(
												  PROPOSAL_IDS, ALL_TARGETS, 
												  FILE_TYPE)
		if redownload_data is False:
			existing_roots = _find_existing_files()
			rows_remove = [i for i, row in enumerate(query_products) if \
						   row['obs_id'] in existing_roots]
			print('Excluding ' + str(len(rows_remove)) +  ' files that have' + \
				' already been downloaded.')				
			query_products.remove_rows(rows_remove)
			
		# Download data to DATA_DIR
		if len(query_products) > 0:
			print('{} new files to download'.format(len(query_products)))
			get_wfc3_data_astroquery.download_products(query_products, 
													   output_dir=DATA_DIR)
			# Remove scans, grisms, full frames from downloaded data
			# Can't do this with astroquery so unfortunatley you must 
			# remove files after. Should eventually filter w/ QL database.
			_remove_unwanted_files()
		else:
			print('No new data to download.')

	###########################################################################		
	############################## sort new data ##############################
	###########################################################################	
	if sort_new_data:
		new_data_files = glob.glob(DATA_DIR+'*{}.fits'.format(FILE_TYPE))

		if len(new_data_files) > 0:
			sort_data.sort_data_targname_filt_propid(DATA_DIR, DATA_DIR, 
													 FILE_TYPE, 
												 	 targname_mappings=\
												 	 TARGNAME_MAPPINGS)
		else:
			print('No new data to sort.')

	dat_dirs = []
	for t in ALL_TARGETS:		
		dat_dirs += glob.glob(DATA_DIR+'{}/*'.format(t))
	###########################################################################		
	############################### run LACosmic ##############################
	###########################################################################	
	if cr_reject:
		pass

	###########################################################################		
	############################ aperture photometry ##########################
	###########################################################################	
	
	if run_ap_phot:
		# create log file for files that fail source detection
		with open(OUTPUT_DIR + 'failed_source_detection.dat', 'w') as fl:
			print('Creating ' + OUTPUT_DIR + 'failed_source_detection.dat')
			fl.write('file,targ,filt\n')
		with open(OUTPUT_DIR + 'flux_outliers.dat', 'w') as fl:
			fl.write('rootname, targ, filt, xcenter, ycenter, flux, percent_dif\n')

		for dirr in dat_dirs:
			targ, filt = dirr.split('/')[-2], dirr.split('/')[-1]
			print('Running photometry on {}, {}'.format(targ, filt))
			phot_file_path = OUTPUT_DIR+'photcats/'+'{}_{}.dat'.format(targ, filt)

			# glob for input input files in proposal subdirectories
			files = glob.glob(dirr+'/*/*{}.fits'.format(FILE_TYPE))

			print('{} files to process'.format(len(files)))
			# run source detection / aperture photometry on each file
			phot_tab_all = Table()
			j=1

			for i, f in enumerate(files): 

				print('\n' + '-'*50 + '{} / {}'.format(i, len(files)) + '-'*50)

				hdu = fits.open(f)
				dat, hdr, hdr0 = hdu[1].data, hdu[1].header, hdu[0].header

				################################################################
				####################### Source detection #######################
				################################################################

				# begin by trying to detect single source with a segmentation. 
				# these parameters work well for most of the standard star data

				coo_tab = aperture_phot.detect_sources_segmap(dat, 
												  threshold=30.,
												  npixels=100, 
												  kernel_fwhm = 1.8, 
												  bkgrnd_threshold=False, 
									  			  show=False) 
				source_detect_flag = 0 # sources that were found on first pass

				# if no sources detected, lower threshold 
				if coo_tab == 0:
					source_detect_flag = 1 # sources not detected initially
					print('No sources detected in {}. lowering threshold.'.format(f))
					coo_tab = aperture_phot.detect_sources_segmap(dat, 
												  threshold=15.,
												  npixels=75, 
												  kernel_fwhm = 1.8, 
												  bkgrnd_threshold=False, 
									  			  show=False) 
					# if second pass fails, give up.
					if coo_tab == 0: 
						print('No sources detected in {}. Moving on.'.format(f))
						with open(OUTPUT_DIR+'failed_source_detection.dat', 'a') as ff:
							ff.write(os.path.basename(f)+','+targ+','+filt+'\n')
						if i == 0:
							j=0
						continue

				# if multiple sources detected, find source closest to 
				# RA_TARG and DEC_TARG. distance from detected source to 
				# RA_TARG and DEC_TARG is useful to store for all detections.
				# do this outside the if statement then filter source table.
				targname = targ
				if targ == 'GRW70':
					targname='GRW +70 5824'
				if targ == 'P330E':
					targname='GSC 02581-02323'

				# Proper motions not accurate in pre-2015, so RA_TARG and 
				# DEC_TARG are off by more than pointing error. Query SIMBAD
				# and apply PM at date of observation. 
				ra_new, dec_new = utils.apply_proper_motion_targ(targname, 
										hdu[0].header['expstart'])
				# Convert to image coordinates
				wcs = WCS(hdr, hdu)
				ras, decs = wcs.all_pix2world(coo_tab['xcentroid'], 
											 coo_tab['ycentroid'], 1)
				c = SkyCoord(ra=ra_new*u.degree, dec=dec_new*u.degree)
				catalog = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
				idx, d2d, d3d = c.match_to_catalog_sky(catalog)

				if len(coo_tab) > 1:
					print('{} sources detected in {}.'.format(len(coo_tab), f))
					source_detect_flag = 2 # multiple sources detected initially
					print('Closest match to target at distance of {}.'.format(d2d))
					coo_tab = coo_tab[idx:idx+1]

				if pamcorr:
					print('Applying PAM to {}.'.format(f))
					dat = utils.make_PAMcorr_file_UVIS(copy.deepcopy(dat), 
													  hdu[0].header, 
													  hdr, PAM_PATH)

				#convert data to countrate 
				exptime = hdu[0].header['exptime']
				dat = dat/exptime
			
				#calculate sky and rms
				xc, yc = coo_tab['xcentroid'], coo_tab['ycentroid']
				r_in, r_out = SKY_ANNULUS
				back, back_std = aperture_phot.calc_sky_annulus(dat, xc, yc, 
													r_in, r_out, 
													sky_method='mode', 
													sigma_clip=SIGMA_CLIP_SKY)
				#aperture photometry. pass in data - back
				phot_tab = aperture_phot.circular_aperture_photometry(dat-back, 
										coo_tab['xcentroid'], 
										coo_tab['ycentroid'],
										aperture_radii=CIRCULAR_APERTURE_RADII)
				print(phot_tab)

				#add other values to phot table
				phot_tab = _add_hdr_values_to_phottab(phot_tab, hdr0)
				phot_tab['back'] = back
				phot_tab['back_std'] = back_std
				phot_tab['source_detect_flag'] = source_detect_flag
				phot_tab['d2_radec_targ'] = Angle(d2d).degree

				#add phot errors 
				for rad in CIRCULAR_APERTURE_RADII:
					flux = phot_tab['countrate_{}'.format(rad)]*exptime
					phot_ap_area = np.pi * rad**2.
					sky_ap_area = (np.pi * r_out**2.) - (np.pi * r_in**2.)
					flux_error = aperture_phot.compute_phot_err_daophot(flux, 
														  back, back_std, 
														  phot_ap_area,
														  sky_ap_area, 
														  gain=1.0)
					phot_tab['err_{}'.format(rad)] = flux_error
				if (i == 0) or (j == 0):
					phot_tab_all = phot_tab
					j = 1
				else:
					phot_tab_all = vstack([phot_tab_all, phot_tab])


			phot_tab = _clip_and_log_flux_outliers(phot_tab_all, 5, targ, filt)

			#write table to file
			print('Writing {}'.format(phot_file_path))
			phot_tab_all.write(phot_file_path, format='csv',overwrite=True)

if __name__ == '__main__':

	main_standard_star_phot_pipeline()