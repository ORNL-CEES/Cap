#!/usr/bin/env bash

set -x
cd ${PREFIX}/source
git clone https://github.com/dalg24/cap-data.git
apt-get install -y clang-format-3.7
cp /usr/bin/clang-format-3.7 /usr/bin/clang-format
apt-get remove -y \
    libpython-stdlib \
    libpython2.7-minimal \
    libpython2.7-stdlib \
    python \
    python-minimal \
    python2.7 \
    python2.7-minimal
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
