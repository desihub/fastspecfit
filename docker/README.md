Build a Docker container for fastspecfit.
=========================================

Using podman-hpc at NERSC
-------------------------

First, build the container:
```
time podman-hpc build --tag fastspecfit:3.1.1 --file ./Containerfile ./
```
then 'migrate' the image so it can be used for jobs:
```
podman-hpc migrate fastspecfit:3.1.1
```

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
podman-hpc run --rm -it fastspecfit:3.1.1 /bin/bash
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
