Getting started
===============

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

cap has a required dependency on C++11.

+------------------------+------------+---------+
| Packages               | Dependency | Version |
+========================+============+=========+
| Boost                  | Required   | 1.59.0  |
+------------------------+------------+---------+
| deal.II                | Optional   | 8.3.0   |
+------------------------+------------+---------+
| GNU Scientific Library | Optional   | 1.16    |
+------------------------+------------+---------+
| Python                 | Optional   | 2.7     |
+------------------------+------------+---------+

boost
^^^^^
Boost version 1.59.0 or later is required.
Boost can be downloaded from `here <http://www.boost.org/users/download>`_.
Assuming that you have downloaded `boost_1_59_0.tar.bz2` into the
`${PREFIX}/archive` directory, boost may be installed by running:

.. code::

    $ mkdir ${PREFIX}/source/boost
    $ tar -xf ${PREFIX}/archive/boost_1_59_0.tar.bz2 -C ${PREFIX}/source/boost --strip-components=1
    $ cd ${PREFIX}/source/boost && ./bootstrap.sh --prefix=${PREFIX}/install/boost
    $ ./b2 install -j<N> variant=release cxxflags="-std=c++11"

deal.II
^^^^^^^
The open source finite element library deal.II is optional.
It is only required to work with energy storage devices of type ``SuperCapacitor``.
Version 8.3.0 or later compiled with C++11 support is required.
The development sources can be found `here <https://github.com/dealii/dealii>`_.

To download the release version 8.3.0, do:

.. code::

    $ wget -output-document=${PREFIX}/archive/dealii-8.3.0.tar.gz https://github.com/dealii/dealii/releases/download/v8.3.0/dealii-8.3.0.tar.gz
    $ mkdir ${PREFIX}/source/dealii && tar -xf ${PREFIX}/archive/dealii-8.3.0.tar.gz -C ${PREFIX}/source/dealii --strip-components=1

It is a good idea to make a `configure_dealii` script such as:

.. code-block:: bash
    :linenos:

    EXTRA_ARGS=$@

    cmake                                                \
        -D CMAKE_INSTALL_PREFIX=${PREFIX}/install/dealii \
        -D CMAKE_BUILD_TYPE=Release                      \
        -D DEAL_II_WITH_CXX11=ON                         \
        -D BOOST_DIR=${PREFIX}/install/boost             \
        $EXTRA_ARGS                                      \ 
        ${PREFIX}/source/dealii

Then run:

.. code::

    $ cd ${PREFIX}/build/dealii
    $ ../configure_dealii
    $ make -j<N> install

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

    $ git clone https://github.com/dalg24/cap.git ${PREFIX}/source/cap
    $ git clone https://github.com/dalg24/cap-data.git ${PREFIX}/source/cap-data

`cap-data` contains a series of 2-D and 3-D meshes to model batteries or supercapacitors.

Create a `configure_cap` script in `${PREFIX}/build`:

.. code-block:: bash
    :linenos:

    EXTRA_ARGS=$@

    cmake                                               \
        -D CMAKE_INSTALL_PREFIX=${PREFIX}/install/cap   \
        -D BOOST_INSTALL_DIR=${PREFIX}/install/boost    \
        -D DEAL_II_INSTALL_DIR=${PREFIX}/install/dealii \
        -D CAP_DATA_DIR=${PREFIX}/source/cap-data       \
        $EXTRA_ARGS                                     \ 
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


Enable the python wrappers
--------------------------

To build the Python wrappers cap must be configured with an extra flag
``PYTHON_INSTALL_DIR`` that tells cmake where Python is installed.

Find out where Python is installed:

.. code::

    $ export PYTHON_INSTALL_DIR=`python -c "import sys; print sys.prefix"`

Configure cap to build the python interface and (re)install:

.. code::

    $ cmake -DPYTHON_INSTALL_DIR=${PYTHON_INSTALL_DIR} ${PREFIX}/source/cap


Build this documentation
------------------------

Run the configuration script with the extra flag:

.. code::

    $ ../configure_cap -DENABLE_DOCUMENTATION=ON

Open the file `index.html` in the directory `docs/html`.
