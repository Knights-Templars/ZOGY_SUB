# Import the necessary packages

import os
import numpy
from astropy.io import fits
import glob
from astropy.coordinates import SkyCoord
from astropy import units as u

# Object Details

OBJECT_NAME = 'SN2017hpa'
RA_hms = '04:39:50.734'
DEC_hms = '+07:03:55.22'
HOST_GALAXY = 'UGC 3122'

def ra_dec(ra, dec):
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    RA = c.ra.degree
    DEC = c.dec.degree
    return RA,DEC


ra, dec = ra_dec(RA_hms, DEC_hms)
RA = round(ra,5)
DEC = round(dec, 5)
print RA, DEC


cwd = os.getcwd()

files = glob.glob('*.fits')
print files

def remove_file(file_name):
    try:
        os.remove(file_name)
    except OSError:
        pass
        
# function for removing files having similar names

def remove_similar_files(common_text):
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)

def group_similar_files(text_list, common_text, exceptions=''):
    list_files=glob.glob(common_text)
    if exceptions !='':
		list_files=filter(lambda z: not re.search(exceptions, z), list_files)

    list_files.sort()
    if len(text_list) !=0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name+'\n')
                
    return list_files  

def text_list_to_python_list(text_list):
    
    if os.path.exists(text_list):
        with open(text_list, 'r+') as f:
            python_list = f.read.split()
            
    return python_list


def python_list_to_text_list(python_list, text_list):
    
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element)+'\n')


#for text in ['*.list']:
    #remove_similar_files(text)
    
def edit_header(files_list):
    
    for file_name in files_list:
        image = fits.open(file_name, 'update')
        file_header = image[0].header
        list_keywords = ['OBJECT', 'RA', 'DEC']
        dict_header = {'OBJECT': OBJECT_NAME, 'RA': RA, 'DEC': DEC}
        comment = {'OBJECT': 'The Name of the Object', 'RA': 'The Right Ascension', 'DEC': 'The Declination'}
    
        for keyword in list_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all = True)
            file_header.append(card=(keyword, dict_header[keyword], comment[keyword]))
            
        image.flush()
        image.close()
        
        
edit_header(files)


           
