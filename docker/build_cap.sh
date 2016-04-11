#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=2}
# fetch the repo with the meshes
git clone https://github.com/dalg24/cap-data.git ${PREFIX}/source/cap-data
# append the option flag --allow-run-as-root to mpiexec
mv ${MPI_DIR}/bin/mpiexec ${MPI_DIR}/bin/mpiexec.alias
echo '#!/usr/bin/env bash\nmpiexec.alias --allow-run-as-root "$@"' > ${MPI_DIR}/bin/mpiexec
chmod +x ${MPI_DIR}/bin/mpiexec
# setup report to cdash
if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
    DASHBOARD=Continuous
else
    DASHBOARD=Experimental
fi
export LD_LIBRARY_PATH=${BOOST_DIR}/lib:${LD_LIBRARY_PATH}
echo 'set(CTEST_DROP_METHOD "http")' >>  ${PREFIX}/source/cap/CTestConfig.cmake
echo 'set(CTEST_DROP_SITE "jupyterdocker.ornl.gov/CDash")' >> ${PREFIX}/source/cap/CTestConfig.cmake
echo 'set(CTEST_DROP_LOCATION "/submit.php?project=Cap")' >> ${PREFIX}/source/cap/CTestConfig.cmake
echo 'set(CTEST_DROP_SITE_CDASH TRUE)' >>  ${PREFIX}/source/cap/CTestConfig.cmake
cat ${PREFIX}/source/cap/CTestConfig.cmake
# build the code
mkdir ${PREFIX}/build/cap
cd ${PREFIX}/build/cap
cmake \
    -G "Unix Makefiles" \
    -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_CXX_FLAGS="-Wall -Wextra" \
    -D BUILD_SHARED_LIBS=ON \
    -D ENABLE_PYTHON=ON \
    -D PYTHON_LIBRARY=${PYTHON_DIR}/lib/libpython3.5.so \
    -D PYTHON_INCLUDE_DIR=${PYTHON_DIR}/include/python3.5 \
    -D BOOST_DIR=${BOOST_DIR} \
    -D ENABLE_DEAL_II=ON \
    -D DEAL_II_DIR=${DEAL_II_DIR} \
    -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
    -D ENABLE_COVERAGE=ON \
    -D ENABLE_FORMAT=ON \
    -D SITE="travis-ci" \
    -D BUILDNAME="\"${CTEST_BUILD_NAME}\"" \
    ${PREFIX}/source/cap
make -j${NPROC} -i
# run unit tests
ctest -j${NPROC} --dashboard ${DASHBOARD} -T Test
# check code coverage
make coverage-cpp && make coverage-python
sed -i.fixpath "s|/dummy/cap||g" lcov.info
sed -i.fixpath "s|/dummy/pycap||g" lcov.info
sed -i.fixpath "s|python/pycap|python/source|g" coverage.xml
