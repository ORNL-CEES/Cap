#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=2}
# append the option flag --allow-run-as-root to mpiexec
cat > /usr/local/bin/mpiexec <<\EOF
#!/usr/bin/env bash
/usr/bin/mpiexec --allow-run-as-root "$@"
EOF
chmod +x /usr/local/bin/mpiexec
# make a couple aliases
ln -s /usr/bin/python3.5 /usr/local/bin/python
ln -s /usr/bin/clang-format-3.7 /usr/local/bin/clang-format
# build the code
mkdir ${PREFIX}/build/cap
cd ${PREFIX}/build/cap
${CAP_DIR}/scripts/docker_cmake -D ENABLE_COVERAGE=ON
make -j${NPROC} -i
# run unit tests
export LD_LIBRARY_PATH=${BOOST_DIR}/lib:${LD_LIBRARY_PATH}
ctest -j${NPROC} -V
# check code coverage
make coverage-cpp
make coverage-python
sed -i.fixpath "s|/dummy/cap||g" lcov.info
sed -i.fixpath "s|/dummy/pycap||g" lcov.info
sed -i.fixpath "s|python/pycap|python/source|g" coverage.xml
