# Build a Podman container for FastSpecFit.

## Overview

In production, we use a custom-built [Podman](https://podman.io) container. The
instructions below illustrate how we build the container for use at
[NERSC](https://docs.nersc.gov), which makes use of
[podman-hpc](https://docs.nersc.gov/development/containers/podman-hpc/overview).

All images are tagged and have been publicly checked into
[dockerhub/desihub](https://hub.docker.com/orgs/desihub/repositories).

## Log into dockerhub (optional)

When building a container, first log into `dockerhub` (credentials required):
```
podman-hpc login docker.io
```

## Build the base container

We first build a "base" container to hold our installation of
[MPICH](https://www.mpich.org/) and
[mpi4py](https://mpi4py.readthedocs.io/en/stable/). Using a base image allows us
to make changes to the top-level container without having to rebuild the base
image, which can take up to 20 minutes. (Note that a two-stage build does not
work at NERSC because the cache is not persistent between login nodes.)

Build, tag, migrate, and push the base container to `dockerhub`:
```
podman-hpc build --tag desihub/fastspecfit-base:1.0 --file ./Containerfile-base .
podman-hpc migrate desihub/fastspecfit-base:1.0
podman-hpc push desihub/fastspecfit-base:1.0
```

## Build the production container

The production container pulls in tagged versions of the requisite Python
dependencies and also pre-compiles and caches all the `numba` functions in the
`FastSpecFit` code base. As before, build, tag, migrate, and push the container:
```
podman-hpc build --tag desihub/fastspecfit:3.1.1 --file ./Containerfile .
podman-hpc migrate desihub/fastspecfit:3.1.1
podman-hpc push desihub/fastspecfit:3.1.1
```

## Deploy the container

To deploy the container in production please refer to the [FastSpecFit
documentation](https://fastspecfit.readthedocs.io/en/latest/) for the latest
details and instructions. However, briefly, one can test the container in an
interactive node:
```
salloc -N 1 -C cpu -A desi -t 01:00:00 --qos interactive
podman-hpc pull desihub/fastspecfit:3.1.1
```

First, make sure `mpi4py` works with the [Cray-optimized version of
MPICH](https://docs.nersc.gov/development/containers/podman-hpc/overview/#using-cray-mpich-in-podman-hpc) by running
```
podman-hpc run --rm --mpi desihub/fastspecfit:3.1.1 python -m mpi4py --mpi-lib-version
MPI VERSION    : CRAY MPICH version 8.1.22.12 (ANL base 3.4a2)
MPI BUILD INFO : Wed Nov 09 12:31 2022 (git hash cfc6f82)
```
and
```
srun --ntasks=4 podman-hpc run --rm --mpi desihub/fastspecfit:3.1.1 python -m mpi4py.bench helloworld
Hello, World! I am process 0 of 4 on nid200010.
Hello, World! I am process 1 of 4 on nid200010.
Hello, World! I am process 2 of 4 on nid200010.
Hello, World! I am process 3 of 4 on nid200010.
```

Next, try running a couple healpixels:
```
srun --ntasks=8 podman-hpc run --rm --mpi --group-add keep-groups --volume=/dvs_ro/cfs/cdirs:/dvs_ro/cfs/cdirs \
  --volume=/global/cfs/cdirs:/global/cfs/cdirs --volume=$PSCRATCH:/scratch desihub/fastspecfit:3.1.1 mpi-fastspecfit \
  --specprod=loa --survey=sv3 --program=dark --healpix=26278,26279 --ntargets=16 --mp=4 --outdir-data=/scratch/fasttest
```

## Handy commands

* List the available images:
```
podman-hpc images
REPOSITORY                          TAG         IMAGE ID      CREATED         SIZE        R/O
localhost/desihub/fastspecfit       3.1.1       31a5c5041d1d  50 minutes ago  1.67 GB     true
localhost/desihub/fastspecfit-base  1.0         abd6485d4cb5  16 hours ago    719 MB      true
```

* Check the installed versions of `mpich` and `mpi4py`:
```
podman-hpc run --rm desihub/fastspecfit:3.1.1 python -m mpi4py --version
mpi4py 4.0.1
```
and
```
podman-hpc run --rm desihub/fastspecfit:3.1.1 python -m mpi4py --mpi-lib-version
MPICH Version:3.4.3
MPICH Release date:Thu Dec 16 11:20:57 CST 2021
MPICH ABI:13:12:1
MPICH Device:ch4:ofi
MPICH configure:--with-device=ch4:ofi
MPICH CC:gcc    -O2
MPICH CXX:g++   -O2
MPICH F77:gfortran -fallow-argument-mismatch  -O2
MPICH FC:gfortran -fallow-argument-mismatch  -O2
```

* To make sure the `numba` cache is being used correctly, in the production
example, above, simply set the `NUMBA_DEBUG_CACHE` environment variable
on-the-fly:
```
srun --ntasks=8 podman-hpc run --rm --mpi --group-add keep-groups --volume=/dvs_ro/cfs/cdirs:/dvs_ro/cfs/cdirs \
  --volume=/global/cfs/cdirs:/global/cfs/cdirs --volume=$PSCRATCH:/scratch --env NUMBA_DEBUG_CACHE=1 desihub/fastspecfit:3.1.1 mpi-fastspecfit \
  --specprod=loa --survey=sv3 --program=dark --healpix=26278,26279 --ntargets=16 --mp=4 --outdir-data=/scratch/fasttest
```

* To delete a container:
```
podman-hpc rmi desihub/fastspecfit:3.1.1
```
