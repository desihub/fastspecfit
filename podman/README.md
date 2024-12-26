# Build a Podman container for FastSpecFit.

## Overview

In production, we use a custom-built [Podman](https://podman.io) container. The
instructions below illustrate how we build the container for use at
[NERSC](https://docs.nersc.gov), which makes use of
[podman-hpc](https://docs.nersc.gov/development/containers/podman-hpc/overview).

All images are tagged and have been publicly checked into
[dockerhub/desihub](https://hub.docker.com/orgs/desihub/repositories).

## Logging into dockerhub (optional)

When building a container, first log into `dockerhub`:
```
podman-hpc login docker.io
Username:
Password:
```

## Build the base container

We first build a "base" container to hold our installation of
[mpich](https://www.mpich.org/) and
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

## Handy commands

* List the available images:
```
podman-hpc images
```

* Check the installed versions of `mpich` and `mpi4py`:
```
podman-hpc run --rm desihub/fastspecfit:3.1.1 python -m mpi4py --version
mpi4py 4.0.1

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

* To delete a container:
```
podman-hpc rmi desihub/fastspecfit:3.1.1
```

```
srun --ntasks=4 podman-hpc run --rm --mpi desihub/fastspecfit-base:1.0 python -m mpi4py.bench helloworld
podman-hpc run --userns keep-id --group-add keep-groups --rm --volume=/dvs_ro/cfs/cdirs:/dvs_ro/cfs/cdirs --volume=/pscratch/sd/i/ioannis:/scratch desihub/fastspecfit:3.1.1 /bin/bash
podman-hpc run --userns keep-id --group-add keep-groups --rm --volume=/dvs_ro/cfs/cdirs:/dvs_ro/cfs/cdirs --volume=/global/cfs/cdirs:/global/cfs/cdirs --volume=/pscratch/sd/i/ioannis:/scratch --env NUMBA_CACHE_DIR=/scratch/numba_cache desihub/fastspecfit:3.1.1 /bin/bash
```
