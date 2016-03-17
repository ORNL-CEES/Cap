.. _installation:

Installation
============

This section provide guidelines for installing cap and its TPLs.
We recommend out-of-source builds.

.. code::

    $ export PREFIX=/path/to/some/working/directory/for/the/project
    $ mkdir -p ${PREFIX}
    $ mkdir ${PREFIX}/archive
    $ mkdir ${PREFIX}/source
    $ mkdir ${PREFIX}/build
    $ mkdir ${PREFIX}/install


Install third-party libraries
-----------------------------

Cap has a required dependency on C++14.

+-----------------------------+------------+---------+
| Packages                    | Dependency | Version |
+=============================+============+=========+
| MPI                         | Required   |         |
+-----------------------------+------------+---------+
| Python                      | Optional   |         |
+-----------------------------+------------+---------+
| Boost                       | Required   | 1.59.0  |
+-----------------------------+------------+---------+
| deal.II with p4est/Trilinos | Optional   | 8.4.0   |
+-----------------------------+------------+---------+
| GNU Scientific Library      | Optional   |         |
+-----------------------------+------------+---------+

Message Passing Interface (MPI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cap should be working with any of the MPI implementations. It has only been
tested with Open MPI, MPICH, and Intel MPI.

Boost
^^^^^
Boost version 1.59.0 or later is required.
Boost can be downloaded from `here <http://www.boost.org/users/download>`_.
Make sure to install **all** the libraries.
Do not forget to add the ``using mpi ;`` directive to the
`project-config.jam` file before building boost.

Assuming that you have downloaded `boost_1_59_0.tar.bz2` into the
`${PREFIX}/archive` directory, boost may be installed by running:

.. code::

    $ mkdir ${PREFIX}/source/boost
    $ tar -xf ${PREFIX}/archive/boost_1_59_0.tar.bz2 -C ${PREFIX}/source/boost --strip-components=1
    $ cd ${PREFIX}/source/boost && ./bootstrap.sh --prefix=${PREFIX}/install/boost
    $ echo "using mpi ;" >> project-config.jam
    $ ./b2 install -j<N> variant=release cxxflags="-std=c++14"

deal.II
^^^^^^^

The open source finite element library deal.II is optional.
It is only required to work with energy storage devices of type ``SuperCapacitor``.
Version 8.4.0 or later compiled with C++14/MPI/Boost/p4est/Trilinos support is required.
The development sources can be found `here <https://github.com/dealii/dealii>`_.

To download the release version 8.4.0, do:

.. code::

    $ wget --output-document=${PREFIX}/archive/dealii-8.4.0.tar.gz \
        https://github.com/dealii/dealii/releases/download/v8.4.0/dealii-8.4.0.tar.gz
    $ mkdir ${PREFIX}/source/dealii && tar -xf ${PREFIX}/archive/dealii-8.4.0.tar.gz \
        -C ${PREFIX}/source/dealii --strip-components=1

It is a good idea to make a `configure_dealii` script such as:

.. code-block:: bash
    :linenos:

    EXTRA_ARGS=$@

    cmake                                                \
        -D CMAKE_INSTALL_PREFIX=${PREFIX}/install/dealii \
        -D CMAKE_BUILD_TYPE=Release                      \
        -D DEAL_II_WITH_CXX14=ON                         \
        -D DEAL_II_WITH_MPI=ON                           \
        -D BOOST_DIR=${PREFIX}/install/boost             \
        -D P4EST_DIR=${PREFIX}/install/p4est             \
        -D TRILINOS_DIR=${PREFIX}/install/trilinos       \
        $EXTRA_ARGS                                      \ 
        ${PREFIX}/source/dealii

Then run:

.. code::

    $ mkdir ${PREFIX}/build/dealii && cd ${PREFIX}/build/dealii
    $ ../configure_dealii
    $ make -j<N> install

Please refer to the deal.II documentation to see how to install
`p4est <https://dealii.org/developer/external-libs/p4est.html>`_ and
`Trilinos <https://dealii.org/developer/external-libs/trilinos.html>`_.

gsl
^^^
The GNU Scientific Library (GSL) is optional.
It is only required for electrochemical impedance spectroscopy to perform
fourier analysis. Note that you don't actually need it if you plan on using
the Python wrappers since the FFT algorithms from NumPy can be leveraged in
that case.

.. code::

    $ mkdir ${PREFIX}/build/gsl
    $ cd ${PREFIX}/build/gsl
    $ ${PREFIX}/source/gsl/configure --prefix=${PREFIX}/install/gsl
    $ make -j<N> install


Install cap from source
-----------------------
Get the source:

.. code::

    $ git clone https://github.com/ORNL-CEES/cap.git ${PREFIX}/source/cap
    $ git clone https://github.com/dalg24/cap-data.git ${PREFIX}/source/cap-data

`cap-data` contains a series of 2-D and 3-D meshes to model batteries or supercapacitors.

Create a `configure_cap` script in `${PREFIX}/build`:

.. code-block:: bash
    :linenos:

    EXTRA_ARGS=$@

    cmake \
        -D CMAKE_INSTALL_PREFIX=${PREFIX}/install/cap \
        -D BOOST_DIR=${PREFIX}/install/boost \
        -D ENABLE_DEAL_II=ON \
        -D DEAL_II_DIR=${PREFIX}/install/dealii \
        -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
        $EXTRA_ARGS \ 
        ${PREFIX}/source/cap

Configure, build and install:

.. code::

    $ mkdir ${PREFIX}/build/cap
    $ cd ${PREFIX}/build/cap
    $ ../configure_cap
    $ make -j<N> && make install


Run the tests:

.. code::

    $ ctest -j<N>


Enable the Python wrappers
--------------------------

To build the Python wrappers cap must be configured with an extra flag
``ENABLE_PYTHON=ON``. It is recommended to use Python 3.X but PyCap has
been successfully built with Python 2.X in the past.

.. code::

    $ ../configure_cap -DENABLE_PYTHON=ON
    $ make install

Prepend the `cap/python` directory to the environment variable `PYTHONPATH`
in order to import the pycap module from your Python interpreter.

.. code::

    $ export PYTHONPATH=${PREFIX}/install/cap/lib/pythonX.Y/site-packages:${PYTHONPATH}

``X.Y`` stands for the version of Python that was used to build PyCap, 
for example 2.7 or 3.5.

Launch Python and try:

.. testcode::

    >>> import pycap
    >>> help(pycap)

A number of Python packages are required to use pycap. We recommend you use
pip to install them:

.. code::

    $ pip install numpy scipy matplotlib cython h5py mpi4py


Build this documentation
------------------------

Run the configuration script with the extra flag:

.. code::

    $ ../configure_cap -DENABLE_DOCUMENTATION=ON

Open the file `index.html` in the directory `docs/html`.
