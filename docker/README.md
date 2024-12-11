Build a Docker container for fastspecfit.
=========================================

Using podman-hpc at NERSC
-------------------------

podman-hpc pull ubuntu:22.04

First, build the container:
```
podman-hpc build -t desihub/fastspecfit:latest .
```

podman-hpc migrate desihub/fastspecfit:latest


podman-hpc login docker.io
podman-hpc push docker.io/desihub/fastspecfit:latest
podman-hpc pull docker.io/desihub/fastspecfit:latest

podman-hpc run --rm -it localhost/desihub/fastspecfit:latest /bin/bash



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

docker buildx build --platform linux/amd64,linux/arm64/v8 --push -t desihub/fastspecfit:3.1.2 .
docker buildx build --platform linux/amd64,linux/arm64/v8 --push -t desihub/fastspecfit:latest .
```

To enter the container (with a shell prompt) on a laptop do:
```
docker pull desihub/fastspecfit:latest
docker run -it desihub/fastspecfit:latest
```
or
```
docker pull desihub/fastspecfit:3.1.2
docker run -it desihub/fastspecfit:3.1.2
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
