FROM ubuntu:16.04

ARG NPROC=4
ARG BUILD_TYPE=Release

RUN apt-get update && apt-get install -y \
        build-essential \
        gfortran \
        wget \
        git \
        cmake \
        lcov \
        valgrind \
        cppcheck \
        clang-format-3.7 \
        libssl-dev \
        libpng-dev \
        libfreetype6-dev \
        libxft-dev \
        libsqlite3-dev \
        libbz2-dev \
        libatlas-base-dev \
        zlib1g-dev \
        libopenmpi-dev \
        libhdf5-dev \
        python3.5-dev \
        python3-tk \
        pandoc \
        npm \
        nodejs \
        nodejs-legacy \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set up time zone
RUN echo "America/New_York" > /etc/timezone && \
    dpkg-reconfigure --frontend noninteractive tzdata
    
ENV PREFIX=/scratch
RUN mkdir -p ${PREFIX} && \
    cd ${PREFIX} && \
    mkdir archive && \
    mkdir source && \
    mkdir build && \
    mkdir install

# install boost
RUN export BOOST_URL=http://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.bz2 && \
    export BOOST_SHA256=beae2529f759f6b3bf3f4969a19c2e9d6f0c503edcb2de4a61d1428519fcb3b0 && \
    export BOOST_ARCHIVE=${PREFIX}/archive/boost_1_63_0.tar.bz2 && \
    export BOOST_SOURCE_DIR=${PREFIX}/source/boost/1.63.0 && \
    export BOOST_BUILD_DIR=${PREFIX}/build/boost/1.63.0 && \
    export BOOST_INSTALL_DIR=/opt/boost/1.63.0 && \
    wget --quiet ${BOOST_URL} --output-document=${BOOST_ARCHIVE} && \
    echo "${BOOST_SHA256} ${BOOST_ARCHIVE}" | sha256sum -c && \
    mkdir -p ${BOOST_SOURCE_DIR} && \
    tar -xf ${BOOST_ARCHIVE} -C ${BOOST_SOURCE_DIR} --strip-components=1 && \
    cd ${BOOST_SOURCE_DIR} && \
    ./bootstrap.sh \
        --prefix=${BOOST_INSTALL_DIR} \
        --with-python=python3.5 \
        && \
    echo "using mpi ;" >> project-config.jam && \
    ./b2 install -j${NPROC} \
        link=shared \
        variant=release \
        --build-dir=${BOOST_BUILD_DIR} \
        && \
    rm -rf ${BOOST_ARCHIVE} && \
    rm -rf ${BOOST_BUILD_DIR} && \
    rm -rf ${BOOST_SOURCE_DIR}

ENV BOOST_DIR=/opt/boost/1.63.0
#ENV LD_LIBRARY_PATH=${BOOST_DIR}/lib:${LD_LIBRARY_PATH}

# install TRILINOS
RUN export TRILINOS_VERSION=12.10.1 && \
    export TRILINOS_URL=https://github.com/trilinos/Trilinos/archive/trilinos-release-12-10-1.tar.gz && \
    export TRILINOS_ARCHIVE=${PREFIX}/archive/trilinos-${TRILINOS_VERSION}.tar.xz && \
    export TRILINOS_SOURCE_DIR=${PREFIX}/source/trilinos/${TRILINOS_VERSION} && \
    export TRILINOS_BUILD_DIR=${PREFIX}/build/trilinos/${TRILINOS_VERSION} && \
    export TRILINOS_INSTALL_DIR=/opt/trilinos/${TRILINOS_VERSION} && \
    wget --quiet ${TRILINOS_URL} --output-document=${TRILINOS_ARCHIVE} && \
    mkdir -p ${TRILINOS_SOURCE_DIR} && \
    tar -xf ${TRILINOS_ARCHIVE} -C ${TRILINOS_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${TRILINOS_BUILD_DIR} && \
    cd ${TRILINOS_BUILD_DIR} && \
    cmake \
        -D CMAKE_INSTALL_PREFIX:PATH=${TRILINOS_INSTALL_DIR} \
        -D CMAKE_BUILD_TYPE:SRTING=RELEASE \
        -D BUILD_SHARED_LIBS:BOOL=ON \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
        -D TPL_ENABLE_MPI:BOOL=ON \
        -D MPI_BASE_DIR:PATH=${MPI_DIR} \
        -D TPL_ENABLE_BLAS:BOOL=ON \
        -D TPL_ENABLE_LAPACK:BOOL=ON \
        -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
        -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
        -D Trilinos_ENABLE_Amesos:BOOL=ON \
        -D Trilinos_ENABLE_AztecOO:BOOL=ON \
        -D Trilinos_ENABLE_Epetra:BOOL=ON \
        -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
        -D Trilinos_ENABLE_Ifpack:BOOL=ON \
        -D Trilinos_ENABLE_ML:BOOL=ON \
        -D Trilinos_ENABLE_MueLu:BOOL=ON \
        -D Trilinos_ENABLE_Sacado:BOOL=ON \
        -D Trilinos_ENABLE_Teuchos:BOOL=ON \
        ${TRILINOS_SOURCE_DIR} && \
    make -j${NPROC} install && \
    rm -rf ${TRILINOS_ARCHIVE} && \
    rm -rf ${TRILINOS_BUILD_DIR} && \
    rm -rf ${TRILINOS_SOURCE_DIR}

ENV TRILINOS_DIR=/opt/trilinos/12.10.1

# install p4est
RUN export P4EST_VERSION=1.1 &&\
    export P4EST_URL=http://p4est.github.io/release/p4est-${P4EST_VERSION}.tar.gz && \
    export P4EST_ARCHIVE=${PREFIX}/archive/p4est-${P4EST_VERSION}.tar.gz && \
    export P4EST_SOURCE_DIR=${PREFIX}/source/p4est/${P4EST_VERSION} && \
    export P4EST_BUILD_DIR=${PREFIX}/build/p4est/${P4EST_VERSION} && \
    export P4EST_INSTALL_DIR=/opt/p4est/${P4EST_VERSION} && \
    mkdir -p ${P4EST_SOURCE_DIR} && \
    wget --quiet ${P4EST_URL} --output-document=${P4EST_ARCHIVE} && \
    tar -xf ${P4EST_ARCHIVE} -C ${P4EST_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${P4EST_BUILD_DIR} && \
    cd ${P4EST_BUILD_DIR} && \
    ${P4EST_SOURCE_DIR}/configure \
        --prefix=${P4EST_INSTALL_DIR} \
        --enable-mpi \
        --enable-shared \
        --disable-static \
        --without-blas \
        CC=mpicc \
        CFLAGS="-O2 -fPIC" \
        CXX=mpicxx \
        F77=mpif77 \
        FC=mpifort \
    && \
    make -j${NPROC} install && \
    rm -rf ${P4EST_ARCHIVE} && \ 
    rm -rf ${P4EST_BUILD_DIR} && \
    rm -rf ${P4EST_SOURCE_DIR}

ENV P4EST_DIR=/opt/p4est/1.1

# install dealii
RUN export DEAL_II_VERSION=8.5.1 && \
    export DEAL_II_URL=https://github.com/dealii/dealii/releases/download/v${DEAL_II_VERSION}/dealii-${DEAL_II_VERSION}.tar.gz && \
    export DEAL_II_SHA1=fe8e4748c1d9650487fb5145429a58e8509de276 && \
    export DEAL_II_ARCHIVE=${PREFIX}/archive/dealii-${DEAL_II_VERSION}.tar.gz && \
    export DEAL_II_SOURCE_DIR=${PREFIX}/source/dealii/${DEAL_II_VERSION} && \
    export DEAL_II_BUILD_DIR=${PREFIX}/build/dealii/${DEAL_II_VERSION} && \
    export DEAL_II_INSTALL_DIR=/opt/dealii/${DEAL_II_VERSION} && \
    wget --quiet ${DEAL_II_URL} --output-document=${DEAL_II_ARCHIVE} && \
    echo "${DEAL_II_SHA1} ${DEAL_II_ARCHIVE}" | sha1sum -c && \
    mkdir -p ${DEAL_II_SOURCE_DIR} && \
    tar -xf ${DEAL_II_ARCHIVE} -C ${DEAL_II_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${DEAL_II_BUILD_DIR} && \
    cd ${DEAL_II_BUILD_DIR} && \
    cmake \
        -D CMAKE_INSTALL_PREFIX=${DEAL_II_INSTALL_DIR} \
        -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -D DEAL_II_WITH_MPI=ON \
        -D DEAL_II_WITH_CXX11=ON \
        -D CMAKE_Fortran_COMPILER=mpifort \
        -D CMAKE_CXX_COMPILER=mpicxx \
        -D CMAKE_C_COMPILER=mpicc \
        -D BOOST_DIR=${BOOST_DIR} \
        -D P4EST_DIR=${P4EST_DIR} \
        -D TRILINOS_DIR=${TRILINOS_DIR} \
        ${DEAL_II_SOURCE_DIR} && \
    make -j${NPROC} install && \
    rm -rf ${DEAL_II_ARCHIVE} && \
    rm -rf ${DEAL_II_BUILD_DIR} && \
    rm -rf ${DEAL_II_SOURCE_DIR}

ENV DEAL_II_DIR=/opt/dealii/8.5.1

ENV LD_LIBRARY_PATH=${LAPACK_DIR}/lib:${LD_LIBRARY_PATH}

# install python requirements using pip
RUN wget https://bootstrap.pypa.io/get-pip.py && \
    python3.5 get-pip.py && \
    rm get-pip.py && \
    pip3.5 install -U pip --no-cache-dir && \
    pip3.5 install --no-cache-dir \
        numpy \
        scipy \
        matplotlib \
        cython \
        pypandoc \
        jupyter \
        ipyparallel \
        ipywidgets \
        mpi4py \
        h5py \
        sympy \
        coverage \
        autopep8 \
        docopt \
        jupyterhub \
        && \
    npm install -g configurable-http-proxy && rm -rf ~/.npm

# install tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.10.0/tini && \
    echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - && \
    mv tini /usr/local/bin/tini && \
    chmod +x /usr/local/bin/tini
