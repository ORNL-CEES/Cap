#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# append the option flag --allow-run-as-root to mpiexec
mv ${MPI_DIR}/bin/mpiexec ${MPI_DIR}/bin/mpiexec.alias
echo '#!/usr/bin/env bash' > ${MPI_DIR}/bin/mpiexc
echo 'mpiexec.alias --allow-run-as-root "$@"' >> ${MPI_DIR}/bin/mpiexec
chmod +x ${MPI_DIR}/bin/mpiexec
# build the code
rm -rf build
mkdir build
cd build
#mkdir ${PREFIX}/build/cap
#cd ${PREFIX}/build/cap
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic -Weffc++"\
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D BUILD_SHARED_LIBS=ON \
    -D ENABLE_PYTHON=ON \
    -D PYTHON_LIBRARY=${PYTHON_DIR}/lib/libpython3.5.so \
    -D PYTHON_INCLUDE_DIR=${PYTHON_DIR}/include/python3.5 \
    -D BOOST_DIR=${BOOST_DIR} \
    -D ENABLE_DEAL_II=ON \
    -D DEAL_II_DIR=${DEAL_II_DIR} \
    -D ENABLE_COVERAGE=OFF \
    -D ENABLE_FORMAT=ON \
    ..
#    ${PREFIX}/source/cap
make -j${NPROC} -i
# run unit tests
export LD_LIBRARY_PATH=${BOOST_DIR}/lib:${LD_LIBRARY_PATH}
ctest -j${NPROC} --no-compress-output -T Test
