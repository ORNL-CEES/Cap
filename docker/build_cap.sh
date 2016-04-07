#!/usr/bin/env bash

set -x
cd ${PREFIX}/source
git clone https://github.com/dalg24/cap-data.git
useradd -m -s /bin/bash -N -u 1000 jovyan
chown jovyan ${PREFIX} -R
su jovyan <<EOF
set -x
export LD_LIBRARY_PATH=${BOOST_DIR}/lib:$LD_LIBRARY_PATH
export PATH=${PYTHON_DIR}/bin:${PATH}
ctest -j2 -V -S ${PREFIX}/source/cap/docker/TravisCI.cmake
cd ${PREFIX}/build/cap
make coverage-cpp && make coverage-python
sed -i.fixpath "s|/dummy/cap||g" lcov.info
sed -i.fixpath "s|python/pycap|python/source|g" coverage.xml
EOF
