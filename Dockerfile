FROM dalg24/cap-stack

# install cap
RUN cd ${PREFIX}/source && \
    git clone https://github.com/dalg24/cap.git && \
    git clone https://github.com/dalg24/cap-data.git && \
    mkdir -p ${PREFIX}/build/cap && \
    cd ${PREFIX}/build/cap && \
    cmake \
        -D CMAKE_INSTALL_PREFIX=/opt/cap \
        -D CMAKE_Fortran_COMPILER=mpifort \
        -D CMAKE_CXX_COMPILER=mpicxx \
        -D CMAKE_C_COMPILER=mpicc \
        -D CMAKE_CXX_FLAGS="-fPIC" \
        -D BOOST_INSTALL_DIR=/opt/boost/1.59.0 \
        -D DEAL_II_INSTALL_DIR=/opt/dealii/8.3.0 \
        -D PYTHON_INSTALL_DIR=/opt/python/2.7.10 \
        -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
        ${PREFIX}/source/cap && \
   make install

ENV PYTHONPATH=/opt/cap/python:${PYTHONPATH}
