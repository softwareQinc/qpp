import subprocess
import sys

from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

p = subprocess.Popen(
    "cmake -S pyqpp -B pyqpp/build",
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)

prefix_eigen = "Detected Eigen3 in: "
prefix_pybind11 = "Detected pybind11 in: "
eigen_path = None
pybind11_path = None

print("Running cmake -S pyqpp -B pyqpp/build")
cmake_output = p.stdout.read().decode("ascii").split("\n")
for line in cmake_output:
    print(line)
    pos_eigen = line.find(prefix_eigen)
    pos_pybind11 = line.find(prefix_pybind11)
    if pos_eigen != -1:
        start = pos_eigen + len(prefix_eigen)
        end = line.find(" ", start)
        if end == -1:  # no space found, take the rest of the line
            eigen_path = line[start:]
        else:
            eigen_path = line[start:end]
    if pos_pybind11 != -1:
        pybind11_path = line[pos_pybind11 + len(prefix_pybind11) :]

if eigen_path is None:
    raise Exception("Eigen3 not found!")

if pybind11_path is None:
    raise Exception("pybind11 not found!")

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
