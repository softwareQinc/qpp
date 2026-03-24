# Base image with minimal Python 3.11 environment
FROM python:3.11-slim-bookworm

# Install system dependencies required for building Quantum++ and Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libeigen3-dev \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for safer execution
RUN useradd -m -c "softwareQ" sq

# Set working directory
WORKDIR /home/sq

# Copy Python dependencies
COPY --chown=sq:sq requirements.txt ./

# Switch to non-root user
USER sq

# Create virtual environment and install Python dependencies + pyqpp
RUN python3 -m venv /home/sq/venv && \
    /home/sq/venv/bin/pip install --no-cache-dir --upgrade pip && \
    /home/sq/venv/bin/pip install --no-cache-dir -r requirements.txt && \
    /home/sq/venv/bin/pip install --no-cache-dir git+https://github.com/softwareqinc/qpp

# Clone Quantum++ source and build qpp_qasm example
RUN git clone --depth 1 --branch main https://github.com/softwareqinc/qpp /home/sq/qpp && \
    cmake -S /home/sq/qpp -B /home/sq/qpp/build && \
    cmake --build /home/sq/qpp/build --target qpp_qasm

# Install Quantum++ system-wide and expose qpp_qasm binary
USER root
RUN cmake --install /home/sq/qpp/build && \
    cp /home/sq/qpp/build/examples/qpp_qasm /usr/local/bin/qpp_qasm

# Clean build artifacts and prepare Jupyter workspace
USER sq
RUN rm -rf /home/sq/qpp/build && mkdir -p /home/sq/notebooks
COPY --chown=sq:sq notebooks/ /home/sq/notebooks/

# Ensure virtual environment is used by default
ENV PATH="/home/sq/venv/bin:$PATH"

# Expose Jupyter port and launch notebook server
EXPOSE 8888
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--notebook-dir=/home/sq/notebooks"]
