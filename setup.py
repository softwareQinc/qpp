import subprocess
import sys
from setuptools import setup
from libs.pybind11.setup_helpers import Pybind11Extension

p = subprocess.Popen("cmake pyqpp",
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)

prefix = "Detecting Eigen3 - done (in "
eigen_path = None
print("Running cmake pyqpp")
for line in p.stdout.read().decode('ascii').split('\n'):
    print(line)
    pos = line.find(prefix)
    if pos != -1:
        eigen_path = line[pos + len(prefix):-1]
        break

if eigen_path is None:
    raise Exception('Eigen3 not found!')

ext_modules = [
    Pybind11Extension(
        "pyqpp",
        ["pyqpp/qpp_wrapper.cpp"],
        extra_compile_args=["-Ilibs", "-Iinclude", "-Iqasmtools/include",
                            "-isystem" + eigen_path],
        cxx_std=17,
        include_pybind11=False,
    ),
]

setup(
    name='pyqpp',
    version='4.0.1',
    description='Python 3 wrapper for Quantum++',
    long_description=open('pyqpp/README.md').read(),
    long_description_content_type='text/markdown',
    author='softwareQ',
    author_email='info@softwareq.ca',
    url='https://github.com/softwareQinc/qpp',
    license='MIT',
    platforms=sys.platform,
    install_requires=[
        'numpy',
    ],
    ext_modules=ext_modules)
