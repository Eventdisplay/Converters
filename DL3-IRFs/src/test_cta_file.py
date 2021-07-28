# Test program for OGADF scheme using ogadf_schema
"""
   usage:  python src/test_cta_file.py <file.fits.gz>
"""
from ogadf_schema.irfs import AEFF_2D, EDISP_2D, RAD_MAX
from ogadf_schema.irfs import PSF_TABLE, PSF_3GAUSS, PSF_KING
from ogadf_schema.irfs import BKG_3D, BKG_2D
from astropy.io import fits
import logging
import sys

if len(sys.argv) != 2:
    raise ValueError('please provide FITS file name for testing')

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s| %(message)s'))
logging.getLogger().addHandler(handler)

f = fits.open(sys.argv[1])
print(*[repr(hdu.header['EXTNAME']) for hdu in f[1:]])

print('Checking effective area HDU')
AEFF_2D.validate_hdu(f['EFFECTIVE AREA'], onerror='log')

print('Checking energy dispersion HDU')
EDISP_2D.validate_hdu(f['ENERGY DISPERSION'], onerror='log')

print('Checking gamma-ray point-spread function HDU')
PSF_3GAUSS.validate_hdu(f['POINT SPREAD FUNCTION'], onerror='log')

print('Checking background HDU')
BKG_3D.validate_hdu(f['BACKGROUND'], onerror='log')
