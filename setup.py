from setuptools import setup
from libs.pybind11.setup_helpers import Pybind11Extension
import subprocess

p = subprocess.Popen("cmake pyqpp",
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)

prefix = "Detecting Eigen3 - done in "
eigen_path = None
for line in p.stdout.read().decode('ascii').split('\n'):
    if line.startswith(prefix):
        eigen_path = line[len(prefix):]
        break

if eigen_path is None:
    raise Exception('Eigen3 not found!')

ext_modules = [
    Pybind11Extension(
        "pyqpp",
        ["pyqpp/qpp_wrapper.cpp"],
        extra_compile_args=["-Ilibs", "-Iinclude", "-Iqasmtools/include",
                            "-I" + eigen_path],
        cxx_std=17,
        include_pybind11=False,
    ),
]

setup(
    name='pyqpp',
    version='3.1',
    description='Python wrapper for Quantum++',
    author='softwareQ',
    author_email='info@softwareq.ca',
    url='https://github.com/softwareQinc/qpp',
    install_requires=[
        'numpy',
    ],
    ext_modules=ext_modules)
