# Test program for OGADF scheme using ogadf_schema
"""
   usage:  python src/test_cta_file.py <file.fits.gz>
           (only 3D implemented at this point)
"""
import logging
import sys

from astropy.io import fits
from ogadf_schema.irfs import AEFF_2D, BKG_3D, EDISP_2D, PSF_3GAUSS

if len(sys.argv) != 2:
    raise ValueError("please provide FITS file name for testing")

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(levelname)s| %(message)s"))
logging.getLogger().addHandler(handler)

print("Testing", sys.argv[1])

f = fits.open(sys.argv[1])
print(*[repr(hdu.header["EXTNAME"]) for hdu in f[1:]])

print("Checking effective area HDU")
AEFF_2D.validate_hdu(f["EFFECTIVE AREA"], onerror="log")

print("Checking energy dispersion HDU")
EDISP_2D.validate_hdu(f["ENERGY DISPERSION"], onerror="log")

print("Checking gamma-ray point-spread function HDU")
PSF_3GAUSS.validate_hdu(f["POINT SPREAD FUNCTION"], onerror="log")

print("Checking background HDU")
BKG_3D.validate_hdu(f["BACKGROUND"], onerror="log")
