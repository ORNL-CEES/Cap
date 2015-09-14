Getting started
===============

Get the source
--------------

.. code::

    $ export PREFIX=/path/to/some/working/directory/for/the/project
    $ mkdir -p ${PREFIX}/source


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

.. code::

    $ ./bootstrap.sh --prefix=${PREFIX}/install/boost
    $ ./b2 install -j<N> variant=release cxxflags="-std=c++11"

deal.II
^^^^^^^
The open source finite element library deal.II is optional.
It is only required to work with energy storage devices of type
``SuperCapacitor``.
Version 8.3.0 or later compiled with C++11 support is required.
The development sources can be found
`here <https://github.com/dealii/dealii>`_.

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
        ${PREFIX}/source/dealii-8.3.0

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

[Recommended] Create a `configure_cap` script:

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

Build this documentation
------------------------

Run the configuration script with the extra flag:

.. code::

    $ ../configure_cap -DENABLE_DOCUMENTATION=ON

Open the file `index.html` in the directory `docs/html`.
