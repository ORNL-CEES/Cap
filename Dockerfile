FROM dalg24/cap-stack

# install cap and run the tests
RUN cd ${PREFIX}/source && \
    git clone https://github.com/ORNL-CEES/Cap.git cap && \
    git clone https://github.com/dalg24/cap-data.git && \
    mkdir -p ${PREFIX}/build/cap && \
    cd ${PREFIX}/build/cap && \
    cmake \
        -D CMAKE_INSTALL_PREFIX=/opt/cap \
        -D BUILD_SHARED_LIB=ON \
        -D MPI_INSTALL_DIR=${MPI_DIR} \
        -D PYTHON_INSTALL_DIR=${PYTHON_DIR} \
        -D BOOST_INSTALL_DIR=${BOOST_DIR} \
        -D DEAL_II_INSTALL_DIR=${DEAL_II_DIR} \
        -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
        ${PREFIX}/source/cap && \
   make -j2 install && \
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

COPY docker/start-notebook.sh /usr/local/bin/
COPY docker/jupyter_notebook_config.py /home/$NB_USER/.jupyter/
RUN chown -R $NB_USER:users /home/$NB_USER/.jupyter
