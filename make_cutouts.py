#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Created on Wed Dec 30 2019
    
    @author: hk
    """

from astropy.io import fits
import glob
import os
import numpy as np
import numpy.ma as ma
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import argparse
from astropy.nddata.utils import Cutout2D
import subprocess
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.time import Time
import argparse
from astropy.nddata.utils import Cutout2D
from photutils import Background2D, MedianBackground
import image_registration
import scipy
from PIL import Image
import requests
from io import BytesIO


def make_cutout(image, size_x, size_y, x_pixels, y_pixels, i):
    hdu = fits.open(image)[0]
    rawData = hdu.data
    rawHeader = hdu.header
    wcs = WCS(rawHeader)
    for x in x_pixels:
        for y in y_pixels:
            print(x,y)
            cutout = Cutout2D(rawData, position=(x, y), size = (size_x, size_y), wcs = wcs)
            hdu.data = cutout.data
            hdu.header.update(cutout.wcs.to_header())
            keyword_to_remove=('BP_0_0', 'BP_0_1','BP_0_2', 'A_3_0', 'B_3_0','BP_3_0', 'B_1_2', 'B_1_1',
                              'B_2_1', 'B_2_0', 'A_ORDER', 'B_0_3', 'B_0_2', 'BP_0_3','B_ORDER', 'BP_ORDER',
                              'BP_1_2', 'AP_ORDER', 'AP_3_0', 'A_1_1', 'BP_2_0', 'A_1_2', 'AP_2_1', 'AP_2_0',
                              'A_0_2', 'A_0_3', 'BP_1_1', 'BP_1_0', 'A_2_0', 'A_2_1', 'AP_1_0','AP_1_1', 'AP_1_2',
                              'BP_2_1', 'AP_0_1', 'AP_0_0', 'AP_0_3', 'AP_0_2')
            for key in keyword_to_remove:
                try:
                    hdu.header.remove(key)
                except:
                    continue
            new_file_name='cut_'+str(i)+'.wcs.proc.fits'
            print(new_file_name)
            newfile=image.replace('wcs.proc.cr.fits', new_file_name)
            hdu.writeto(newfile, overwrite=True)
            i+=1
    return i


def get_cutout(image, size_x, size_y):
#    x_pixels=[513, 1540, 2567, 3594]
#    y_pixels=[513, 1540, 2567, 3594]
    x_pixels=[700, 2100, 3550]
    y_pixels=[700, 2100, 3550]
    i=1
    i=make_cutout(image, size_x, size_y, x_pixels, y_pixels, i)
    #i=make_cutout(image, 1800, 1800, [2050], [2050], i)
    

def get_filelist(path):
    path=os.getcwd()
    os.chdir(path)
    filelist = glob.glob('*proc.cr.fits')
    print(filelist)
    return filelist
    
    
    
def get_panstarrs_cutouts(filename):
    hdu = fits.open(filename)[0]
    Data = hdu.data
    Header = hdu.header
    filt=Header['FILTER']
    [x,y]=(1500,1500)
    print([x,y])
    wcs = WCS(Header)
    [ra, dec]= wcs.all_pix2world([x],[y],0)
    ra = float('%.5f'%(ra[0]))
    dec=float('%.5f'%(dec[0]))
    print(ra,dec)
    Header['cent_ra']=ra
    Header['cent_dec']=dec
    hdu.header=Header
    hdu.writeto(filename, overwrite=True)
    #[x_img, y_img] = wcs.all_world2pix(ra, dec, 1)
    #x = int(x_img)
    #y=int(y_img)
#         if "_17" in filename:
#             imsize=25
#         else:
    imsize=25
    ps1_command = ('panstamps', '--width=%.3f'%imsize, '--filters=%s'%filt, "--downloadFolder=./panstarrs/", 'stack', '%.5f'%ra, '%.5f'%dec)
    print(ps1_command)
    try:
        rval = subprocess.call(ps1_command)
        print(rval)
    except subprocess.CalledProcessError as err:
        print('Error running panstamps. Exit')
        sys.exit(1)
    
            
def get_legacy_cutouts(filename, ref_list):
    hdu = fits.open(filename)[0]
    Data = hdu.data
    Header = hdu.header
    filt=Header['FILTER']
    [x,y]=(750,750)
    print([x,y])
    wcs = WCS(Header)
    [ra, dec]= wcs.all_pix2world([x],[y],0)
    ra = float('%.5f'%(ra[0]))
    dec=float('%.5f'%(dec[0]))
    print(ra,dec)
    legacy_file='legacy_'+filt+'_'+str(ra)+'_'+str(dec)+'.fits'
    Header['cent_ra']=ra
    Header['cent_dec']=dec
    if legacy_file in ref_list:
        print("Reference file already exists. Skipping file download")
        return None
    else:
        hdu.header=Header
        hdu.writeto(filename, overwrite=True)
        imsize=2000
        url='http://legacysurvey.org/viewer/fits-cutout?ra={}&dec={}&size={}&layer=dr8&pixscale=0.67&bands={}'.format(ra,dec,imsize,filt)
        print(url)
        req = requests.get(url)
        image = Image.open(BytesIO(req.content))
        data = fits.open(url)
        
        print(legacy_file)
    #    legacy_path=path+'/legacy/'
        data.writeto(legacy_file, overwrite=True)



def rename_pan_images(path, filename):
    pan_path=path+'/panstarrs/'
    os.chdir(pan_path)
    panstarrs_list=glob.glob('*stack*.fits')
    os.chdir(default_path)
    hdu = fits.open(filename)[0]
    header = hdu.header
    ra=header['cent_ra']
    dec=header['cent_dec']
    print(ra,dec)
    for panfile in panstarrs_list:
        if str(ra) in panfile:
            print(ra)
            print(filename)
            new_panfilename= 'pan_'+ filename
            print('file found {}. Renaming file to {} '.format(panfile, new_panfilename))
            os.rename(pan_path+panfile, pan_path+new_panfilename)
            print(rval)



    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument("-ra",help="RA in decimal degrees")
#    parser.add_argument("-dec",help="Dec in decimal degrees.")
    parser.add_argument("-path",help="path to the files")
    args = parser.parse_args()
    path=args.path
    
    size_x=1500; size_y=1500
#    out_file='/home/growth/sub_cut_phot_table.csv'
#    f=open(out_file,'w+')
#    f.write('obj,JD,filter_img,mag,mag_err,lim_mag\n')
    os.chdir(path)
    filelist=get_filelist(path)
    print(filelist)

    for file_name in filelist:
        cut_file=file_name.split('-')[0]
        print(cut_file)
        cut_file_phrase="*"+str(cut_file)+"*cut_*.fits"
        print(cut_file_phrase)
        get_cutout(file_name,size_x, size_y)
        cutoutlist=glob.glob(cut_file_phrase)
        print(cutoutlist)
        ref_list=glob.glob('*legacy*.fits')
        print(ref_list)
        for file in cutoutlist:
            #try:
            get_legacy_cutouts(file, ref_list)
            print("Downloading legacy reference image. This step may take some time. Please be patient.")
#            except:
#                get_panstarrs_cutouts(file)
#                print("Downloading panstarrs/legacy reference image. This step may take some time. Please be patient.")
#                rename_pan_images(path, file)
#            science_file=file
#            psf_file=science_file.replace('_resamp.fits','.wcs.proc.fits')
#            psf_file=psf_file.replace('sub_wcs_', 'moffat_')
#            header_file=psf_file.replace('moffat_','')
#            psf_file=psf_file.replace('.fits','.fits.fits')
#            print(psf_file, science_file,header_file)
#            try:
#                filter_img, JD, mag, mag_err, lim_mag =do_psf(psf_file, science_file, header_file)
#                f.write('{},{},{},{},{},{}\n'.format(obj, JD, filter_img, mag, mag_err, lim_mag))
#                print(obj, JD, filter_img, mag, mag_err,  lim_mag)
#            except:
#                print("Error occured while performing subtracted photmetry")
#    f.close()

