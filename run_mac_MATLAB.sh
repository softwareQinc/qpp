#!/bin/sh
# Runs Quantum++ executable with MATLAB support
# The full path to the executable is specified in the command line

# Modify as needed
MATLAB=/Applications/MATLAB_R2017b.app
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$MATLAB/bin/maci64

# Run the executable 
exec $1
