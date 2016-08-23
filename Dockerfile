FROM dalg24/cap-stack

# TODO: move to pre-built image
RUN  ln -sf python3.5 /usr/bin/python

# install cap and run the tests
RUN cd ${PREFIX}/source && \
    git clone https://github.com/ORNL-CEES/Cap.git cap && \
    mkdir -p ${PREFIX}/build/cap && \
    cd ${PREFIX}/build/cap && \
    cmake \
        -D CMAKE_CXX_COMPILER=mpicxx \
        -D CMAKE_INSTALL_PREFIX=/opt/cap \
        -D CMAKE_BUILD_TYPE=Release \
        -D BUILD_SHARED_LIBS=ON \
        -D ENABLE_PYTHON=ON \
        -D BOOST_DIR=${BOOST_DIR} \
        -D ENABLE_DEAL_II=ON \
        -D DEAL_II_DIR=${DEAL_II_DIR} \
        ${PREFIX}/source/cap && \
   make -j2 install && \
   rm -rf ${PREFIX}/build/cap && \
   rm -rf ${PREFIX}/source/cap && \
   rm -rf ${PREFIX}/source/cap-data

ENV PYTHONPATH=/opt/cap/lib/python3.5/site-packages:${PYTHONPATH}
# TODO: this is a tmp fix until adjustments are made to cmake
ENV LD_LIBRARY_PATH=/opt/cap/lib:${LD_LIBRARY_PATH}

ENV SHELL /bin/bash
ENV NB_USER jovyan
ENV NB_UID 1000

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    mkdir /home/$NB_USER/.jupyter && \
    mkdir -p -m 700 /home/$NB_USER/.local/share/jupyter && \
    mkdir -p /home/$NB_USER/.ipython/profile_mpi && \
    mkdir /home/$NB_USER/notebooks && \
    chown -R $NB_USER:users /home/$NB_USER

RUN mkdir /notebooks && chmod -R 777 /notebooks

EXPOSE 8888
WORKDIR /notebooks
VOLUME /notebooks
ENTRYPOINT ["tini", "--"]
CMD ["start-notebook.sh"]

COPY docker/start-notebook.sh /usr/local/bin/
COPY docker/jupyter_notebook_config.py /home/$NB_USER/.jupyter/
COPY docker/jupyter_nbconvert_config.py /home/$NB_USER/.jupyter/
COPY docker/ipengine_config.py /home/$NB_USER/.ipython/profile_mpi/
COPY docker/ipcluster_config.py /home/$NB_USER/.ipython/profile_mpi/
RUN chown -R $NB_USER:users /home/$NB_USER
