import subprocess
import re
import sys

from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

# --------------------------------------------------------------------------- #
# Run CMake and capture its output
# --------------------------------------------------------------------------- #
proc = subprocess.Popen(
    ["cmake", "-S", "pyqpp", "-B", "pyqpp/build"],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
)

assert proc.stdout is not None
output = proc.stdout.read().splitlines()

# --------------------------------------------------------------------------- #
# Regex patterns
# --------------------------------------------------------------------------- #
# Example expected line:
#   Installed/Found Eigen3 version 5.0.1: /path/to/eigen3
re_eigen = re.compile(r"Eigen3 version \d+\.\d+\.\d+:\s*(\S+)")
# Expected line:
#   pybind11: /path/to/pybind11
re_pybind11 = re.compile(r"pybind11:\s*(\S+)")

eigen_path = None
pybind11_path = None

# --------------------------------------------------------------------------- #
# Parse CMake output
# --------------------------------------------------------------------------- #
for line in output:
    m = re_eigen.search(line)
    if m:
        eigen_path = m.group(1)
    m = re_pybind11.search(line)
    if m:
        pybind11_path = m.group(1)

if eigen_path is None:
    raise Exception("Eigen3 not found in CMake output!")

if pybind11_path is None:
    raise Exception("pybind11 not found in CMake output!")

# --------------------------------------------------------------------------- #
# Build pyqpp
# --------------------------------------------------------------------------- #
source_files = ["pyqpp/qpp_wrapper.cpp"]

ext_modules = [
    Pybind11Extension(
        "pyqpp",
        source_files,
        extra_compile_args=[
            "-I" + pybind11_path + "/include",
            "-I" + eigen_path,
            "-Iinclude",
            "-Iqasmtools/include",
            "-Ipyqpp/include",
            "-DQPP_IDX_DEFAULT",
            "-DQPP_BIGINT_DEFAULT",
            "-DQPP_FP_DEFAULT",
        ],
        cxx_std=17,
        include_pybind11=False,
    ),
]

setup(platforms=sys.platform, ext_modules=ext_modules)
