#!/usr/bin/env bash

cd ${PREFIX}/source && \
git clone https://github.com/dalg24/cap-data.git && \
mkdir -p ${PREFIX}/build/cap && \
cd ${PREFIX}/build/cap && \
cmake \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_INSTALL_PREFIX=/opt/cap \
    -D CMAKE_BUILD_TYPE=Debug \
    -D BUILD_SHARED_LIBS=ON \
    -D ENABLE_COVERAGE=ON \
    -D ENABLE_PYTHON=ON \
    -D BOOST_DIR=${BOOST_DIR} \
    -D PYTHON_LIBRARY=${PYTHON_DIR}/lib/libpython3.5.so \
    -D PYTHON_INCLUDE_DIR=${PYTHON_DIR}/include/python3.5 \
    -D ENABLE_DEAL_II=ON \
    -D DEAL_II_DIR=${DEAL_II_DIR} \
    -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
    ${PREFIX}/source/cap && \
make -j2 && \
useradd -m -s /bin/bash -N -u 1000 jovyan && \
chown jovyan ${PREFIX}/build/cap -R && \
su jovyan <<EOF
export LD_LIBRARY_PATH=${BOOST_DIR}/lib:$LD_LIBRARY_PATH && \
export PATH=${PYTHON_DIR}/bin:${PATH} && \
ctest -j2 -V && \
make coverage-cpp && make coverage-python && \
sed -i.fixpath "s|python/pycap|python/source|g" coverage.xml && \
cd ${PREFIX}/source/cap && \
codecov --disable gcov search pycov --file ${PREFIX}/build/cap/lcov.info ${PREFIX}/build/cap/coverage.xml
EOF
