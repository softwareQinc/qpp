# Docker

A self-explanatory minimalistic Docker file is provided in
[`Dockerfile`](docker/Dockerfile).

Build the image by executing

```shell
docker build -t qpp .
```

Run the Jupyter server in a container by executing

```shell
docker run -it -p8888:8888 qpp
```

In case you want to use the Docker container as a development environment,
mount your directory (in this example the current directory) in a Docker
container with

```shell
docker run --rm -it --workdir=/home/sq/hostdir -v ${PWD}:/home/sq/hostdir qpp /bin/bash
```
