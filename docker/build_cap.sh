#!/usr/bin/env bash

cd ${PREFIX}/source && \
git clone https://github.com/dalg24/cap-data.git && \
mkdir -p ${PREFIX}/build/cap && \
cd ${PREFIX}/build/cap && \
cmake \
    -D CMAKE_INSTALL_PREFIX=/opt/cap \
    -D CMAKE_CXX_FLAGS="-fPIC" \
    -D BOOST_INSTALL_DIR=/opt/boost/1.59.0 \
    -D DEAL_II_INSTALL_DIR=/opt/dealii/8.3.0 \
    -D PYTHON_INSTALL_DIR=/opt/python/2.7.11 \
    -D MPI_INSTALL_DIR=/opt/openmpi/1.10.1 \
    -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
    ${PREFIX}/source/cap && \
make install && \
cd python && \
ctest -V
