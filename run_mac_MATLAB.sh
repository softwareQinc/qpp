#!/bin/sh

MATLAB=/Applications/MATLAB_R2017b.app
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$MATLAB/bin/maci64

./build/qpp
