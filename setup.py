import subprocess
import sys

from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup

p = subprocess.Popen(
    "cmake pyqpp", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
)

prefix = "Detecting Eigen3 - done (in "
eigen_path = None
print("Running cmake pyqpp")
for line in p.stdout.read().decode("ascii").split("\n"):
    print(line)
    pos = line.find(prefix)
    if pos != -1:
        eigen_path = line[pos + len(prefix) : -1]
        break

if eigen_path is None:
    raise Exception("Eigen3 not found!")

source_files = ["pyqpp/qpp_wrapper.cpp"]
ext_modules = [
    Pybind11Extension(
        "pyqpp",
        source_files,
        extra_compile_args=[
            "-Ilibs",
            "-Iinclude",
            "-Iqasmtools/include",
            "-Ipyqpp/include",
            "-I" + eigen_path,
            "-DQPP_IDX_DEFAULT",
            "-DQPP_BIGINT_DEFAULT",
            "-DQPP_FP_DEFAULT",
        ],
        cxx_std=17,
        include_pybind11=False,
    ),
]

setup(platforms=sys.platform, ext_modules=ext_modules)
