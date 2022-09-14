#!/bin/bash

# Build script for test and development system crusher at OLCF
# See https://github.com/QMCPACK/qmcpack/pull/4123 for more details on the module file if needed

echo "Loading QMCPACK dependency modules for crusher"
module load cmake/3.22.2
module load cray-fftw
module load openblas/0.3.17-omp
module load rocm
module load cray-hdf5-parallel
module load fftw
module load boost
module load libxml2
