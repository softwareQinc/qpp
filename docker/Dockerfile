FROM ubuntu:latest
ENV TZ="America/Toronto"
RUN apt-get update && \
    apt-get -y update && \
    apt-get install -yq tzdata && \
    ln -fs /usr/share/zoneinfo/America/Toronto /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

RUN apt-get install -y build-essential python3.10 python3-pip python3-dev \
    libeigen3-dev cmake sudo git vim

RUN pip3 -q install pip --upgrade
RUN pip3 install jupyter matplotlib numpy
RUN pip3 install git+https://github.com/softwareqinc/qpp

# Enable a normal user with sudo access
RUN useradd -m -c "softwareQ" sq
RUN echo '%sq ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER sq

# Clone and install Quantum++
WORKDIR /home/sq
RUN git clone --depth 1 --branch main https://github.com/softwareqinc/qpp
WORKDIR /home/sq/qpp
RUN cmake -B build -DWITH_UNIT_TESTS=ON -DWITH_EXAMPLES=ON && \
    sudo cmake --build build --target install
USER sq
WORKDIR /home/sq

# Create a notebook directory for Jupyter
RUN mkdir notebooks
