#--------------------------------------------------------------------------------------------------#


from astropy.io import fits
import sys
import glob
import os
from astroquery.gaia import Gaia
import subprocess
import numpy as np
import numpy.ma as ma
import warnings
from astropy.wcs import WCS
import astropy.units as u
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.coordinates import SkyCoord
from astropy.io import ascii

#--------------------------------------------------------------------------------------------------#
cwd = os.getcwd()
sys.path.append(cwd)
import ldactools as aw
import py_zogy
warnings.simplefilter('ignore')
astrom_scamp = cwd + '/scamp.conf'
astrom_sex = cwd + '/astrom.sex'
astrom_swarp = cwd + '/config.swarp'
psf_config_name = cwd + '/photom.psfex'
sexpsfparam = cwd + '/sexpsf.param'
astrom_param = cwd + '/astrom.param'

#--------------------------------------------------------------------------------------------------#

fits_delete = False

#--------------------------------------------------------------------------------------------------#

def remove_file(file_name):

    try:
        os.remove(file_name)
    except OSError:
        pass     
    
def remove_similar_files(ctext):

    for file_name in glob.glob(ctext):
        remove_file(file_name)
        
        
for text in ['*.cat', '*.resamp', '*.resamp.bkgrms*', '*.scaled', '*.unc', '*.weight', '*.head', '*.ldac', '*.txt', '*.ps', '*.zip', '*.psf', \
            '*.psfmodel','*.weight.fits']:
    remove_similar_files(ctext = text)   
    
    if fits_delete:
        remove_similar_files(ctext = '*.fits')    
    
#--------------------------------------------------------------------------------------------------#

class color:

   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

        
    
    
def display_text(text):

    print('#'+'-'*(10+len(text))+'#')
    print('#'+('-'*5)+ color.BOLD+ color.UNDERLINE+str(text)+ color.END+('-'*5)+'#')
    print('#'+'-'*(10+len(text))+'#')        

#--------------------------------------------------------------------------------------------------#

def get_ps1_template(ra, dec, imsize, filt):

        display_text("Downloading Ps1 template")

	'''
	ps_list = glob.glob('*stack*skycell*fits')
	[os.remove(p) for p in ps_list]
	ini_file_list = glob.glob('stack*fits')
	
	ps1_command = ['panstamps', '--width=%.3f'%imsize, '--filters=%s'%filt, "--downloadFolder=.", 'stack', '%.5f'%ra, '%.5f'%dec]
	print(ps1_command)
	try:
		rval = subprocess.call(ps1_command)
	except subprocess.CalledProcessError as err:
		print('Error running panstamps. Exit')
		sys.exit(1)
	'''
	new_file_list = glob.glob('stack*fits')
	new_stacks = [x for x in new_file_list]
	
	if len(new_stacks) == 0:
		print('panstamps did not produce anything ..')
		return None
	ps1_filename = new_stacks[0]
	
	print('panstamps created new file %s'%ps1_filename)
	return ps1_filename
	
	
def make_gaia_catalog(ra, dec, catalog_box_size, catalog_min_mag, catalog_max_mag, catname, writetext = False, writeldac = False):

        display_text('Generate Gaia Catalog')
    
	job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.panstarrs1_best_neighbour AS pbest, gaiadr2.panstarrs1_original_valid AS ps1 WHERE g.source_id = pbest.source_id AND pbest.original_ext_source_id = ps1.obj_id AND CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND ps1.r_mean_psf_mag > %.2f AND ps1.r_mean_psf_mag < %.2f AND pmra IS NOT NULL AND pmdec IS NOT NULL AND abs(pmdec) > 0 AND abs(pmdec) < 40 AND abs(pmra)>0 AND abs(pmra) < 40 AND ps1.n_detections > 6 AND pbest.number_of_mates=0 AND pbest.number_of_neighbours=1;"%(ra, dec, catalog_box_size/60, catalog_min_mag, catalog_max_mag), dump_to_file = False)
	
	p = job.get_results()

	p['ra_errdeg'] = p['ra_error'] / 3.6e6
	p['dec_errdeg'] = p['dec_error'] / 3.6e6
	p['FLAGS'] = 0
	
	p.remove_columns(['astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'astrometric_weight_al', 'astrometric_pseudo_colour', 'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations', 'visibility_periods_used', 'astrometric_sigma5d_max', 'frame_rotator_object_type', 'matched_observations', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'l', 'b', 'ecl_lon', 'ecl_lat', 'priam_flags', 'teff_val', 'teff_percentile_lower', 'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper', 'e_bp_min_rp_val', 'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper', 'flame_flags', 'radius_val', 'radius_percentile_lower', 'radius_percentile_upper', 'lum_val', 'lum_percentile_lower', 'lum_percentile_upper', 'gaia_astrometric_params', 'obj_name', 'obj_id', 'ra_2', 'dec_2', 'ra_error_2', 'dec_error_2', 'epoch_mean', 'zone_id', 'obj_info_flag', 'quality_flag', 'designation', 'phot_variable_flag', 'datalink_url', 'epoch_photometry_url', 'original_ext_source_id'])
	
	if writetext:
		ascii.write(p, 'catalog.txt')
	if writeldac:
		if os.path.exists(catname + '.ldac'):
			os.remove(catname + '.ldac')
		aw.save_table_as_ldac(p, catname + '.ldac')
   
		
	return 1
	
def run_sextractor(sciname, refname):

        display_text('Running Sextractor')
	
	sex_command = ['sex', '-c', astrom_sex, sciname, '-CATALOG_NAME', sciname +'.cat', '-WEIGHT_TYPE', 'NONE', '-CHECKIMAGE_TYPE', 'NONE']
	try:
		rval = subprocess.call(sex_command)
	except subprocess.CalledProcessError as err:
		print('Could not run sextractor on %s. Exit'%sciname)
		sys.exit(1)
		
	sex_command = ['sex', '-c', astrom_sex, refname, '-CATALOG_NAME', refname +'.cat', '-WEIGHT_TYPE', 'NONE', '-CHECKIMAGE_TYPE', 'NONE']	
	try:
		rval = subprocess.call(sex_command)
	except subprocess.CalledProcessError as err:
		print('Could not run sextractor on %s. Exit'%refname)
		sys.exit(1)
	
	
def run_scamp(filelistname, ra ,dec, imagesize, sciname, refname):

        display_text('Running Scamp')
	
	catname = filelistname + '.cat'
	make_gaia_catalog(ra, dec, imagesize, 10, 22, catname, writetext = True, writeldac = True)
	
	scamp_command = ['scamp', '-c', astrom_scamp, '@%s'%filelistname, '-ASTREFCAT_NAME', catname + '.ldac']
	
	try:
		rval = subprocess.call(scamp_command)
	except subprocess.CalledProcessError as err:
		print('Could not run scamp on %s. Exit'%refname)
		sys.exit(1)
	
	os.rename(sciname+'.head', sciname[:-5]+'.head')
	os.rename(refname+'.head', refname[:-5]+'.head')
	
	
def run_swarp(scifilename, reffilename, pxscale, imgpixsize, ra, dec):

        display_text('Running Swarp')

	run_sextractor(scifilename, reffilename)
	filelistname = '%s_list.txt'%scifilename
	f = open(filelistname, 'w')
	f.write('%s.cat\n%s.cat'%(scifilename, reffilename))
	f.close()
	
	run_scamp(filelistname, ra, dec, imgpixsize*pxscale/60, scifilename, reffilename)
	
	#Create weight files for science and referece images
	scihdu = fits.open(scifilename, 'update')
	scidata = scihdu[0].data
	refhdu = fits.open(reffilename, 'update')
	refdata = refhdu[0].data
	refhdu.close()
	
	#Create a weight map to throw out the NANs
	sci_weight_data = np.ones(scidata.shape)
	sci_weight_data[np.isnan(scidata)] = 0
	ref_weight_data = np.ones(refdata.shape)
	ref_weight_data[np.isnan(refdata)] = 0
	sci_weight_hdu = fits.PrimaryHDU(sci_weight_data)
	sci_weight_hdu.writeto(scifilename + '.weight', overwrite = True)
	ref_weight_hdu = fits.PrimaryHDU(ref_weight_data)
	ref_weight_hdu.writeto(reffilename + '.weight', overwrite = True)	
	
	#swarp_command = "swarp %s -c %s -RESAMPLE Y -RESAMPLE_DIR . -SUBTRACT_BACK Y -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.3f -IMAGE_SIZE %d -CENTER %.8f,%.8f -COMBINE Y -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s"%(scifilename, astrom_swarp, scifilename + '.weight', pxscale, imgpixsize, ra , dec, scifilename + '.resamp', scifilename + '.resamp.weight')

        swarp_command = "swarp %s -c %s -RESAMPLE Y -RESAMPLE_DIR . -SUBTRACT_BACK Y -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.3f -IMAGE_SIZE %d -CENTER %.8f,%.8f -COMBINE N -WEIGHTOUT_NAME %s"%(scifilename, astrom_swarp, scifilename + '.weight', pxscale, imgpixsize, ra, dec, scifilename+'.resamp.weight')	
	
	print(swarp_command)
	os.system(swarp_command)
	
	#swarp_command = "swarp %s -c %s -RESAMPLE Y -RESAMPLE_DIR . -SUBTRACT_BACK N -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.3f -IMAGE_SIZE %d -CENTER %.8f,%.8f -COMBINE Y -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s"%(reffilename, astrom_swarp, reffilename + '.weight', pxscale, imgpixsize, ra , dec, reffilename + '.resamp', reffilename + '.resamp.weight')
	
	swarp_command = "swarp %s -c %s -RESAMPLE Y -RESAMPLE_DIR . -SUBTRACT_BACK N -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.3f -IMAGE_SIZE %d -CENTER %.8f,%.8f -COMBINE N -WEIGHTOUT_NAME %s"%(reffilename, astrom_swarp, reffilename + '.weight', pxscale, imgpixsize, ra , dec, reffilename+'.resamp.weight')
	
	print(swarp_command)
	os.system(swarp_command)
	
	sciresampname = scifilename[:-5] + '.resamp.fits'
	refresampname = reffilename[:-5] + '.resamp.fits'
	
	scipsf = get_psf(sciresampname)
	#scipsf = 'wcs_20190116_034-g.fits.psfmodel'
	refpsf = get_psf(refresampname)

	scihdu.close()
	
	return sciresampname, refresampname, scipsf, refpsf
	
	
def get_psf(filename):

        display_text('Running Sextractor')

	sex_command = ['sex',  '-c',  astrom_sex, filename, '-CATALOG_NAME', filename+'.cat', '-WEIGHT_TYPE', 'NONE', '-CHECKIMAGE_TYPE', 'BACKGROUND_RMS', '-CHECKIMAGE_NAME', filename+'.bkgrms']
	try:
		rval = subprocess.call(sex_command)
	except subprocess.CalledProcessError as err:
		print('Could not run sextractor on %s. Exit'%filename)
	        sys.exit(1)
	        
        display_text('Ruuning PSFEX')	        
    	
        psfex_command = ['psfex', '-c', psf_config_name, filename+'.cat']
        psfname = filename + '.psf'
	
        try:
    	    rval = subprocess.call(psfex_command)
        except subprocess.CalledProcessError as err:
            print('Could not run psfex on %s. Exit'%filename)
	    sys.exit(1)
	
	psf_model_data = fits.open(psfname)[1].data[0][0][0]
	psf_model_data = psf_model_data / np.sum(psf_model_data)
	psf_model_hdu = fits.PrimaryHDU(psf_model_data)
	psf_model_hdu.writeto(filename + '.psfmodel', overwrite = True)
	
	return psfname+'model'

	
	
def run_zogy(sciname, refname, scipsf, refpsf):

        display_text('Ruuning ZOGY')
	
	#Cross match science and reference image catalogs to get flux scaling factor and astometric uncertainties
	sci_catalog = aw.get_table_from_ldac(sciname + '.cat')
	print sci_catalog
	ref_catalog = aw.get_table_from_ldac(refname + '.cat')
	print ref_catalog
	
	good_sci_sources = (sci_catalog['FLAGS'] == 0)  & (sci_catalog['SNR_WIN'] > 20) & (sci_catalog['FWHM_WORLD'] < 4./3600) & (sci_catalog['FWHM_WORLD'] > 0.5/3600) # & (sci_catalog['SNR_WIN'] < 100)
	good_ref_sources = (ref_catalog['FLAGS'] == 0)  & (ref_catalog['SNR_WIN'] > 20) & (ref_catalog['FWHM_WORLD'] < 4./3600) & (ref_catalog['FWHM_WORLD'] > 0.5/3600) # & (ref_catalog['SNR_WIN'] < 100)
	
	print('Number of good sources SCI: %d REF: %d'%(np.sum(good_sci_sources), np.sum(good_ref_sources)))

	ref_catalog = ref_catalog[good_ref_sources]
	sci_catalog = sci_catalog[good_sci_sources]
	
	#sci_coords = SkyCoord(ra = sci_catalog['ALPHAWIN_J2000'], dec = sci_catalog['DELTAWIN_J2000'], frame = 'icrs')
	sci_coords = SkyCoord(ra = sci_catalog['X_WORLD'], dec = sci_catalog['Y_WORLD'], unit = 'degree')
	ref_coords = SkyCoord(ra = ref_catalog['X_WORLD'], dec = ref_catalog['Y_WORLD'], unit = 'degree')
	
	#Cross match the catalogs
	idx_ref, idx_sci, d2d, d3d = sci_coords.search_around_sky(ref_coords, 1.0*u.arcsec)

	xpos_sci = sci_catalog['XWIN_IMAGE']
	ypos_sci = sci_catalog['YWIN_IMAGE']
	xpos_ref = ref_catalog['XWIN_IMAGE']
	ypos_ref = ref_catalog['YWIN_IMAGE']
	sci_flux_auto = sci_catalog['FLUX_AUTO']
	ref_flux_auto = ref_catalog['FLUX_AUTO']

	print('Number of cross-matched sources is %d'%(len(d2d)))
	
	ast_unc_x = np.std(xpos_sci[idx_sci] - xpos_ref[idx_ref])
	ast_unc_y = np.std(ypos_sci[idx_sci] - ypos_ref[idx_ref])
	
	print('Astrometric uncertainties are X: %.2f Y:%.2f'%(ast_unc_x, ast_unc_y))
	
	flux_scale_mean, flux_scale_median, flux_scale_std = sigma_clipped_stats(sci_flux_auto[idx_sci] / ref_flux_auto[idx_ref])
	flux_scale = flux_scale_median
	print('Flux scaling for reference is %.5f +/- %.5f'%(flux_scale, flux_scale_std))
	
	scihdu = fits.open(sciname)
	sci_header = scihdu[0].header
	sci_data = scihdu[0].data
	sci_gain = 1.6
	scihdu.close()
	sciweighthdu = fits.open(sciname[:-5]+'.weight.fits')
	sci_weight_data = sciweighthdu[0].data
	sciweighthdu.close()
	refhdu = fits.open(refname)
	ref_header =refhdu[0].header
	ref_data = refhdu[0].data
	ref_gain = 1.04
	refhdu.close()
	refweighthdu = fits.open(refname[:-5]+'.weight.fits')
	ref_weight_data = refweighthdu[0].data
	refweighthdu.close()
	
	image_mask = (sci_weight_data == 0) | (ref_weight_data == 0)
	sci_data[image_mask] = 0
	ref_data[image_mask] = 0
	
	ref_data = ref_data * flux_scale
	
	sci_rms = 0.5*(np.percentile(sci_data[sci_data!=0], 84.13) - np.percentile(sci_data[sci_data!=0], 15.86))
	ref_rms = 0.5*(np.percentile(ref_data[ref_data!=0], 84.13) - np.percentile(ref_data[ref_data!=0], 15.86))
	
	print('Science RMS is %.2f. Reference RMS is %.2f'%(sci_rms, ref_rms))

	sci_poisson_noise = np.copy(sci_data) / sci_gain
	sci_poisson_noise[sci_poisson_noise < 0] = 0
	sci_rms_image = np.sqrt(sci_poisson_noise + sci_rms**2)
	ref_poisson_noise = np.copy(ref_data) / ref_gain
	ref_poisson_noise[ref_poisson_noise < 0] = 0
	ref_rms_image = np.sqrt(ref_poisson_noise + ref_rms**2)


	scirmshdu = fits.PrimaryHDU(sci_rms_image)
	scirmshdu.writeto(sciname + '.unc', overwrite=True)
	sci_rms_input = sciname + '.unc'
	refrmshdu = fits.PrimaryHDU(ref_rms_image)
	refrmshdu.writeto(refname + '.unc', overwrite=True)
	ref_rms_input = refname + '.unc'
	sciscaledhdu = fits.PrimaryHDU(sci_data)
	sciscaledhdu.writeto(sciname + '.scaled', overwrite = True)
	sci_diff_input = sciname + '.scaled'
	refscaledhdu = fits.PrimaryHDU(ref_data)
	refscaledhdu.writeto(refname + '.scaled', overwrite = True)
	ref_diff_input = refname + '.scaled'
	
	D, P_D, S_corr = py_zogy.py_zogy(sci_diff_input, ref_diff_input, scipsf, refpsf, sci_rms_input, ref_rms_input, sci_rms, ref_rms, dx=ast_unc_x, dy=ast_unc_y)
	
	diff_hdu = fits.PrimaryHDU(D)
	diff_hdu.header = sci_header
	diff_hdu.writeto(sciname + '.diff', overwrite = True)
	
	diff_psf_hdu = fits.PrimaryHDU(P_D)
	diff_psf_hdu.writeto(sciname + '.diff.psf', overwrite = True)	
	
	scorr_hdu = fits.PrimaryHDU(S_corr)
	scorr_hdu.header = sci_header
	scorr_mean, scorr_median, scorr_std = sigma_clipped_stats(S_corr)
	
	print('Scorr Mean, Median, STD is %.3f, %.3f, %.3f'%(scorr_mean, scorr_median, scorr_std))
	scorr_hdu.writeto(sciname + '.scorr', overwrite = True)
	


def main(scifilename):
	
	imagehdu = fits.open(scifilename, 'update')
	imageHead = imagehdu[0].header
	w = WCS(imageHead)
	cd11 = imageHead['CD1_1']
	cd21 = imageHead['CD2_1']
	cd12 = imageHead['CD1_2']
	cd22 = imageHead['CD2_2']
	nxpix = imageHead['NAXIS1']
	nypix = imageHead['NAXIS2']
	image_x_cent = nxpix/2
	image_y_cent = nypix/2
	
	[raCent, decCent] = w.all_pix2world(image_x_cent, image_y_cent, 1)
	
	xscale = np.sqrt(cd11**2 + cd21**2)
	yscale = np.sqrt(cd12**2 + cd22**2)

	try:	
		filterName = imageHead['FILTER']
		filterName = filterName.strip("'")[0]
	except KeyError:
		print('FILTER keyword not found!')
		filterName = None
	
	boxwidth = max(xscale * nxpix, yscale * nypix)*60*1.2
	
	print('Science image centered at RA %.3f DEC %.3f with SIZE of %.2f arcmin in FILTER %s'%(raCent, decCent, boxwidth, filterName))
	
	#reffilename = get_ps1_template(raCent, decCent, boxwidth, filterName)
	
	#reffilename = 'cutout_157.4405_29.4368.fits'
	
	reffilename = 'pan_' + scifilename
	print reffilename

	sciresampname, refresampname, scipsf, refpsf = run_swarp(scifilename, reffilename, xscale * 3600, max(nxpix, nypix), raCent, decCent)
	
	print sciresampname, refresampname
	
	run_zogy(sciresampname, refresampname, scipsf, refpsf)

	imagehdu.close()
	
	

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description = 'Code to subtract input image from a PS1 template')
	parser.add_argument('imagename', help = 'FITS file to perform subtraction on')
	
	args = parser.parse_args()
	image = args.imagename

	main(image)
	
	
