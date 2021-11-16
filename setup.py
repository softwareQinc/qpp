from glob import glob
from setuptools import setup
from libs.pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "pyqpp",
        ["pyqpp/qpp_wrapper.cpp"],
        extra_compile_args=["-Ilibs", "-Iinclude", "-Iqasmtools/include"],
        cxx_std=17,
        include_pybind11=False,
    ),
]

setup(
    name='pyqpp',
    version='3.0',
    description='Python wrapper for Quantum++',
    author='softwareQ',
    author_email='info@softwareq.ca',
    url='https://github.com/softwareQinc/qpp',
    ext_modules=ext_modules)
