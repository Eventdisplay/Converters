# Eventdisplay to DL2 FITS format

Example of use:

```shell
python generate_DL2_file.py gamma_cone.S.3HB9-FD_ID0.eff-0.root myfile.fits lapalma
```

Note the DL2 format used is based on the open DL3 specifications from:

https://github.com/open-gamma-ray-astro/gamma-astro-data-formats

The only modifications implemented were:
* Adding MC parameters (true energy, true azimuth/zenith...) by adding to the column name (`MC_*`) 
within the events table:

```bash
ColDefs(
    name = 'OBS_ID'; format = '1K'
    name = 'EVENT_ID'; format = '1K'
    name = 'MC_ALT'; format = '1E'; unit = 'deg'
    name = 'MC_AZ'; format = '1E'; unit = 'deg'
    name = 'MC_ENERGY'; format = '1E'; unit = 'TeV'
    name = 'ALT'; format = '1E'; unit = 'deg'
    name = 'AZ'; format = '1E'; unit = 'deg'
    name = 'ENERGY'; format = '1E'; unit = 'TeV'
    name = 'MULTIP'; format = '1J'
)
```

* Adding an additional binary table to the FITS file containing a histogram with the number of simulated events vs `MC_ENERGY`.
```bash
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       7   ()      
  1  EVENTS        1 BinTableHDU     44   549346R x 9C   [1K, 1K, 1E, 1E, 1E, 1E, 1E, 1E, 1J]   
  2  SIMULATED EVENTS    1 BinTableHDU     21   1R x 3C   [60E, 60E, 60E]  

ColDefs(
    name = 'MC_ENERG_LO'; format = '60E'; unit = 'TeV'
    name = 'MC_ENERG_HI'; format = '60E'; unit = 'TeV'
    name = 'EVENTS'; format = '60E'
)
```

 