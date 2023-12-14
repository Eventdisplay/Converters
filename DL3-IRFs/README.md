# IRF conversion from CTAO-ROOT to FITS-GDAF format

see [GDAF website](https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/irf_components/index.html#point-spread-function) for details on the data model.

## Installation

Preferred installation is using docker or apptainer.

Use the prepared docker image to be downloaded from ghcr.io/eventdisplay/converters/dl3-irfs .

Easiest to use it interactively by opening a bash in the corresponding image, e.g.:

```bash
apptainer exec docker://ghcr.io/eventdisplay/converters:latest-dl3-irfs bash
```

Built the docker image on your own:

```bash
docker build -t dl3-irfs .
```

Without docker, all what is needed is a [ROOT](https://root.cern.ch/) and [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) installation, and then compile it with:

```bash
make convertSensitivityFilesToFITS
```

## Converting IRFs in ROOT to FITS

General syntax:

```bash
./convertSensitivityFilesToFITS <file.root> <file.fits.gz> 3D
```

(3D option is the default)

Using the docker image discussed above, e.g.:

```bash
$  docker run --rm -it -v "$(pwd)":/workdir \
   ghcr.io/eventdisplay/converters:latest-dl3-irfs
   /workdir/convertSensitivityFilesToFITS \
   /workdir/Prod5-South-20deg-AverageAz-14MSTs37SSTs.180000s-v0.1.root \
   /workdir/Prod5-South-20deg-AverageAz-14MSTs37SSTs.180000s-v0.1.fits.gz \
   3D
```

## Testing

Install ogadf tools:

```bash
conda env create -f environment.yml
conda activate dl3irfs
```

To test the converted FITS file, do:

```bash
python src/test_cta_file.py <file.fits.gz>
```

## Notebooks

Some testing notebooks are in the notebook directory. Use `jupyter lab` to run them in a gammapy environment.
