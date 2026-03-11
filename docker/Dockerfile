FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libeigen3-dev \
    git \
    vim \
    sudo

# Create user with sudo access in one go
RUN useradd -m -c "softwareQ" sq && \
    echo 'sq ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER sq
WORKDIR /home/sq

COPY --chown=sq:sq requirements.txt ./
RUN python3 -m venv venv && \
    ./venv/bin/pip install --no-cache-dir --upgrade pip && \
    ./venv/bin/pip install --no-cache-dir -r requirements.txt && \
    ./venv/bin/pip install --no-cache-dir git+https://github.com/softwareqinc/qpp

# Clone Quantum++ and build qpp_qasm
RUN git clone --depth 1 --branch main https://github.com/softwareqinc/qpp && \
    cmake -S qpp -B qpp/build && \
    cmake --build qpp/build --target qpp_qasm && \
    sudo cmake --build qpp/build --target install && \
    sudo cp qpp/build/qpp_qasm /usr/local/bin && \
    rm -rf qpp/build  # Clean up build artifacts

# Create a notebook directory for Jupyter
RUN mkdir notebooks
ENV PATH="/home/sq/venv/bin:$PATH"
# CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--notebook-dir=./notebooks"]
