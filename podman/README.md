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
to make changes...


Next, build, migrate, and push the base container:
```
podman-hpc build --tag desihub/fastspecfit-base:1.0 --file ./Containerfile-base ./
podman-hpc migrate desihub/fastspecfit-base:1.0
podman-hpc push desihub/fastspecfit-base:1.0
```

Check the versions:
```
podman-hpc run --rm desihub/fastspecfit-base:1.0 mpirun --version
podman-hpc run --rm desihub/fastspecfit-base:1.0 python -m mpi4py --version
podman-hpc run --rm desihub/fastspecfit-base:1.0 python -m mpi4py --mpi-std-version
podman-hpc run --rm desihub/fastspecfit-base:1.0 python -m mpi4py --mpi-lib-version
```

```
srun --ntasks=4 podman-hpc run --rm --mpi desihub/fastspecfit-base:1.0 python -m mpi4py.bench helloworld
```


### Test scripts:
podman-hpc login docker.io
podman-hpc build --tag desihub/fastspecfit-test --file ./Containerfile-test .
podman-hpc migrate desihub/fastspecfit-test
srun --ntasks=2 podman-hpc run --rm --mpi desihub/fastspecfit-test /usr/local/bin/test-mpi-fastspecfit --mp=4


###

First, build the container:
```
podman-hpc build --tag desihub/fastspecfit:3.1.1 --file ./Containerfile ./
podman-hpc migrate desihub/fastspecfit:3.1.1
podman-hpc push desihub/fastspecfit:3.1.1
```

podman-hpc run --userns keep-id --rm --volume=/dvs_ro/cfs/cdirs:/dvs_ro/cfs/cdirs --volume=/pscratch/sd/i/ioannis:/scratch fastspecfit:3.1.1 /bin/bash

podman-hpc run --userns keep-id --group-add keep-groups --rm --volume=/dvs_ro/cfs/cdirs:/dvs_ro/cfs/cdirs --volume=/global/cfs/cdirs:/global/cfs/cdirs --volume=/pscratch/sd/i/ioannis:/scratch --env NUMBA_CACHE_DIR=/scratch/numba_cache desihub/fastspecfit:3.1.1 /bin/bash

List the available images:
```
podman-hpc images
```


```
podman-hpc login docker.io
podman-hpc push fastspecfit:3.1.1
podman-hpc pull desihub/fastspecfit:3.1.1
```

Enter the container:
```
podman-hpc run --rm fastspecfit:3.1.1 /bin/bash
```

To delete a container:
```
podman-hpc rmi desihub/fastspecfit:3.1.1
```

Legacy Instructions
-------------------

Build a cross-platform docker container as documented [here](https://www.docker.com/blog/faster-multi-platform-builds-dockerfile-cross-compilation-guide), [here](https://blog.jaimyn.dev/how-to-build-multi-architecture-docker-images-on-an-m1-mac/), and [here](https://docs.nersc.gov/development/shifter/how-to-use/).

The very first time, create a builder instance with
```
docker buildx create --name FastSpecFit-build --node FastSpecFit --use
```
and then subsequently simply load that instance with
```
docker buildx use FastSpecFit-build
```

Then, subsequently, to create a new (or the latest) version or tag, do:
```
export DOCKER_BUILDKIT=0
export COMPOSE_DOCKER_CLI_BUILD=0

docker buildx build --platform linux/amd64,linux/arm64/v8 --push -t desihub/fastspecfit:3.1.1 .
docker buildx build --platform linux/amd64,linux/arm64/v8 --push -t desihub/fastspecfit:latest .
```

To enter the container (with a shell prompt) on a laptop do:
```
docker pull desihub/fastspecfit:latest
docker run -it desihub/fastspecfit:latest
```
or
```
docker pull desihub/fastspecfit:3.1.1
docker run -it desihub/fastspecfit:3.1.1
```

Or at NERSC:
```
shifterimg pull docker:desihub/fastspecfit:latest
shifter --image docker:desihub/fastspecfit:latest bash
```

To install the jupyter kernel do:
```
mkdir -p ~/.local/share/jupyter/kernels/fastspecfit
wget -O ~/.local/share/jupyter/kernels/fastspecfit/kernel.json https://raw.githubusercontent.com/desihub/fastspecbi/main/etc/jupyter-kernel.json
```

To grab the setup script (and modify it as needed) do:
```
wget https://raw.githubusercontent.com/desihub/fastspecbi/main/bin/fastspecfit-setup.sh
```
