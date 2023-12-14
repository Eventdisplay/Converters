# Eventdisplay to DL2 FITS format

Example of use
----------

```shell
python ./generate_DL2_file.py -l LAYOUT_NAME gamma_onSource.S.BL-4LSTs25MSTs70SSTs-MSTF_ID0.eff-0.root
```

Format definition
-------------
Note the DL2 format used is based on the open specifications from:

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

DL2 events trees are called 'DL2EventTree' and include gamma/hadron cut statistics.


### Data trees

Data tree variables (subset):

- MCe0: MC energy (in TeV)
- MCaz, MCel: MC shower direction
- MCxoff and MCyoff: true position (relative to camera centre; in deg)
- ArrayPointing_Azimuth, ArrayPointing_Elevation: telescope pointing direction
- erec: reconstructed energy (in TeV)
- nimages: number of images used for the reconstruction
- xoff and xoff: reconstructed position (relative to camera centre; in deg)
- Class: cut class defining different type of events (see below)
- MVA: BDT mva parameter

The Class parameter can be used to select different type of events:

1. Events passing gamma/hadron separation cut and direction cut

```
DL2EventTree->Draw("MVA", "Class==5" );
```

2. Events passing gamma/hadron separation cut and not direction cut

```
DL2EventTree->Draw("MVA", "Class==0" );
```

3. Events before applying gamma/hadron separation cut and before applying direction cut

```
DL2EventTree->Draw("MVA", "Class==0||Class==7||Class==5", "");
```

### Simulated events

The simulated events vs energy are stored in a histogram in the file:

```
TH1D *h = (TH1D*)gDirectory->Get("hEmcUW");
```

Events are not weighted while filling this histogram (in contrary to earlier versions)

Links
-----

CTA internal prod3b redmine page: https://redmine.cta-observatory.org/projects/cta_analysis-and-simulations/wiki/Eventdisplay_Prod3b_DL2_Lists
(older version; needs to be updated)

## Installation

Install the required packages and activate conda environment:

```
mamba env create -f environment.yml
conda activate DL2
```

Update your environment:
```
conda env update -f environment.yml
```

## Docker images

Docker images can be downloaded from the package directory: https://github.com/Eventdisplay/Converters/pkgs/container/converters

Use the container, e.g., by:
```
$ docker run --rm -it -v "$(pwd):/workdir" \
   ghcr.io/eventdisplay/converters:latest-dl2 \
   ./generate_DL2_file.py -l LAYOUT_NAME \
   /workdir/gamma_onSource.S.BL-4LSTs25MSTs70SSTs-MSTF_ID0.eff-0.root
```

Use apptainer:

```
apptainer exec docker://ghcr.io/eventdisplay/converters:latest-dl2 bash
```
