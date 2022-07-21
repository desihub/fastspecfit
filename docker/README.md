Build a Docker container for fastspecfit.
=========================================

Build a cross-platform docker container as documented [here](https://www.docker.com/blog/faster-multi-platform-builds-dockerfile-cross-compilation-guide) and [here](https://blog.jaimyn.dev/how-to-build-multi-architecture-docker-images-on-an-m1-mac/).
```
export DOCKER_BUILDKIT=0
export COMPOSE_DOCKER_CLI_BUILD=0

docker buildx create --use
docker buildx build --platform linux/amd64,linux/arm64/v8 --push -t desihub/fastspecfit:v1.0 .
docker buildx build --platform linux/amd64,linux/arm64/v8 --push -t desihub/fastspecfit:latest .
```

To enter the container (with a shell prompt) on a laptop do:
```
docker pull desihub/fastspecfit:latest
docker run -it desihub/fastspecfit:latest
```
or
```
docker pull desihub/fastspecfit:v1.0
docker run -it desihub/fastspecfit:v1.0
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
