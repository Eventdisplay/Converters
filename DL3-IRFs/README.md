# IRF conversion to GDAF format

see [GDAF website](https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/irf_components/index.html#point-spread-function) for format details.

## Installation

Requires ROOT and cfitsio installation, e.g. by:
```
conda env create -f environment.yml
conda activate dl3irfs
```

Compilation:
```
make convertSensitivityFilesToFITS
```

## Converting IRFs in ROOT to FITS

```
./convertSensitivityFilesToFITS <file.root> <file.fits.gz> 3D
```
(3D option is the default)

To test the converted FITS file, do:
```
python src/test_cta_file.py <file.fits.gz>
```

## Notebooks

Some testing notebooks are in the notebook directory. Use `jupyter lab` to run them in a gammapy environment.
Test IRFs can be found on https://forge.in2p3.fr/projects/cta-s-optimization/wiki

## Docker

Docker images are used only for easier debugging.

### Building

```
$ docker build -t dl3-irfs .
```

### Running

```
$ docker run --rm -it -v "$(pwd):/workdir" dl3-irfs bash
```

### Debugging

Install valgrind:
```
apt-get install g++ valgrind -y
```

Run with suppressions:
```
valgrind --suppressions=/data/root/etc/valgrind-root.supp ./convertSensitivityFilesToFITS
```

