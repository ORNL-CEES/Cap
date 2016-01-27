FROM dalg24/cap-stack

ARG CAP_BRANCH=master
ARG CAP_REPOSITORY=dalg24/cap

# install cap and run the tests
RUN cd ${PREFIX}/source && \
    git clone --branch ${CAP_BRANCH} https://github.com/${CAP_REPOSITORY}.git cap && \
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
   ctest -V && \
   rm -rf ${PREFIX}/build/cap && \
   rm -rf ${PREFIX}/source/cap && \
   rm -rf ${PREFIX}/source/cap-data

ENV PYTHONPATH=/opt/cap/python:${PYTHONPATH}

ENV SHELL /bin/bash
ENV NB_USER jovyan
ENV NB_UID 1000

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    mkdir /home/$NB_USER/.jupyter && \
    mkdir /home/$NB_USER/.local && \
    mkdir /home/$NB_USER/notebooks && \
    chown -R $NB_USER:users /home/$NB_USER

EXPOSE 8888
WORKDIR /home/$NB_USER/notebooks
VOLUME /home/$NB_USER/notebooks
ENTRYPOINT ["tini", "--"]
CMD ["start-notebook.sh"]

COPY start-notebook.sh /usr/local/bin/
COPY jupyter_notebook_config.py /home/$NB_USER/.jupyter/
RUN chown -R $NB_USER:users /home/$NB_USER/.jupyter
