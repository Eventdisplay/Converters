# IRF conversion to GDAF format

see [GDAF website](https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/irf_components/index.html#point-spread-function) for format details.

**this is for now a playground for Gernot to learn more about the format; don't expect lot's of progress and don't think this is for production**

## Converting IRFs in ROOT to FITS

```
./convertSensitivityFilesToFITS
```

## Notebooks

Some testing notebooks are in the notebook directory. Use `jupyter lab` to run them.

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

