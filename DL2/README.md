# Eventdisplay to DL2 FITS format

Example of use
----------

```shell
python ./generate_DL2_file.py gamma_onSource.S.BL-4LSTs25MSTs70SSTs-MSTF_ID0.eff-0.root paranal
```

Format definition
-------------
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

* Adding an additional HDU containing a histogram with the number of simulated events vs `MC_ENERGY`.
```bash
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU       7   ()      
  1  EVENTS        1 BinTableHDU     44   549346R x 9C   [1K, 1K, 1E, 1E, 1E, 1E, 1E, 1E, 1J]   
  2  SIMULATED EVENTS    1 BinTableHDU     21   1R x 3C   [60E, 60E, 60E]  
```
The table is stored as a binary table with 3 columns: two columns defining the low and high `MC_ENERGY` bins, and a third column containing 
the number of events.
```bash
ColDefs(
    name = 'MC_ENERG_LO'; format = '60E'; unit = 'TeV'
    name = 'MC_ENERG_HI'; format = '60E'; unit = 'TeV'
    name = 'EVENTS'; format = '60E'
)
```

Input format from Eventdisplay Effective Area files
----------

## Event trees

DL2 events trees are called 'data' and gamma/hadron cut statistics are listed in 'fEventTreeCuts' (both trees have the same number of events, so the TTree::AddFriend() mechanism can be used). Both trees contain entries for all trees after reconstruction.


### Data trees

Data tree variables (subset):

- MCe0 - MC energy (in TeV)
- MCxcore and MCycore - MC core position (in m)
- MCxoff and MCyoff - true position (relative to camera centre; in deg)
- ErecS - reconstructed energy (in TeV)
- Xcore, YCore - reconstructed core position (in m)
- Xoff and Yoff - reconstructed position (relative to camera centre; in deg)
- NImages - number of images used for the reconstruction

### Gamma/hadron cut trees:

The tree fEventTreeCuts can be used to select different type of events:

1. Events passing gamma/hadron separation cut and direction cut

```
fEventTreeCuts->Draw("MVA", "Class==5" );
```

2. Events passing gamma/hadron separation cut and not direction cut

```
fEventTreeCuts->Draw("MVA", "Class==0" );
```

3. Events before applying gamma/hadron separation cut and before applying direction cut

```
fEventTreeCuts->Draw("MVA", "Class==0||Class==7||Class==5", "");
```

### Simulated events¶

The simulated events vs energy are stored in a histogram in the file:

```
TH1D *h = (TH1D*)gDirectory->Get("hEmcUW");
```

Events are not weighted while filling this histogram (in contrary to earlier versions)

 
Links
-----

CTA internal prod3b redmine page: https://forge.in2p3.fr/projects/cta_analysis-and-simulations/wiki/Eventdisplay_Prod3b_DL2_Lists
