# IRF conversion to GDAF format

see [GDAF website](https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/irf_components/index.html#point-spread-function) for format details.

## Installation

Preferred installation is using docker:

```
$ docker build -t dl3-irfs .
```

Otherwise all what is needed is a ROOT and cfitsio installation, and then compile it with:
```
make convertSensitivityFilesToFITS
```

## Converting IRFs in ROOT to FITS

General syntax:

```
./convertSensitivityFilesToFITS <file.root> <file.fits.gz> 3D
```
(3D option is the default)

Using the docker image discussed above, e.g.:

```
$  docker run --rm -it -v "$(pwd):/workdir" dl3-irfs \
   /data/Converters/DL3-IRFs/convertSensitivityFilesToFITS \
   /workdir/data/Prod5-South-20deg-AverageAz-14MSTs37SSTs.180000s-v0.1.root \
   /workdir/Prod5-South-20deg-AverageAz-14MSTs37SSTs.180000s-v0.1.fits.gz \
   3D
``` 

## Testing

Install ogadf tools:

```
conda env create -f environment.yml
conda activate dl3irfs
```

To test the converted FITS file, do:
```
python src/test_cta_file.py <file.fits.gz>
```

## Notebooks

Some testing notebooks are in the notebook directory. Use `jupyter lab` to run them in a gammapy environment.
Test IRFs can be found on https://forge.in2p3.fr/projects/cta-s-optimization/wiki

