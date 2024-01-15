#!/bin/bash

# Build the C code using Makefile
make ascot5_main MPI=1
make libascot MPI=1

# Copy the compiled C library to the conda environment
cp build/libascot.so $PREFIX/lib/
cp build/ascot5_main $PREFIX/bin/

# Install the Python package
$PYTHON -m pip install . --ignore-installed --prefix=$PREFIX
#python setup.py install ### --no-deps
