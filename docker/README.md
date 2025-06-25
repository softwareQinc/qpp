## Docker

A self-explanatory minimalistic Docker file is provided
in [`Dockerfile`](https://github.com/softwareqinc/qpp/tree/main/docker/Dockerfile).

Build the image by executing

```shell
docker build -t softwareq-qpp .
```

---

Run the Jupyter server in a container by executing

```shell
docker run -p8888:8888 -it --workdir=/home/sq/notebooks softwareq-qpp sh -c ". ~/venv/bin/activate && jupyter notebook --port=8888 --no-browser --ip=0.0.0.0"
```

---

In case you want to use the Docker container as a development environment,
mount your directory (in this example the current directory) in a Docker
container with

```shell
docker run --rm -it --workdir=/home/sq/hostdir -v ${PWD}:/home/sq/hostdir softwareq-qpp /bin/bash
```
