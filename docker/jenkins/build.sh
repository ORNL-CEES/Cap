#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# append the option flag --allow-run-as-root to mpiexec
cat > /usr/local/bin/mpiexec <<\EOF
#!/usr/bin/env bash
#/usr/bin/mpiexec --allow-run-as-root "$@"
EOF
chmod +x /usr/local/bin/mpiexec
# make a couple aliases
ln -s /usr/bin/python3.5 /usr/local/bin/python
ln -s /usr/bin/clang-format-3.7 /usr/local/bin/clang-format
# build the code
cd ${PREFIX}/source/cap
rm -rf build
mkdir build
cd build
#mkdir ${PREFIX}/build/cap
#cd ${PREFIX}/build/cap
cmake \
    -G "Unix Makefiles" \
    -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic -Weffc++"\
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D BUILD_SHARED_LIBS=ON \
    -D ENABLE_PYTHON=ON \
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

#cppcheck \
#    --std=c++11 \
#    --enable=all \
#    --inconclusive \
#    --xml --xml-version=2 \
#    -I ../cpp/source/dummy \
#    -I ../cpp/source/deal.II/dummy \
#    ../cpp 2> cppcheck.xml
