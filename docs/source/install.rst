.. _installation:

Installation
============

This section provide guidelines for installing Cap from source.

Note that it is **not** necessary to build Cap from source to use it. Refer to
the :ref:`Docker<docker>` section for instructions on how to pull the latest
image of Cap.

Third-party libraries
---------------------

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

Cap and its dependencies may be built using `spack
<https://github.com/llnl/spack>`_. You would need to install the following
packages:

.. code::

    $ spack install boost +graph +icu_support +mpi +python
    $ spack install trilinos ~hypre ~mumps +boost \
         ^boost+graph+icu_support+mpi+python
    $ spack install dealii~arpack~gsl~oce~petsc+trilinos+mpi \
         ^trilinos~hypre~mumps ^boost+graph+icu_support+mpi+python
    $ spack install py-mpi4py
    $ spack install py-matplotlib
    $ spack install py-h5py

Before buiding Cap, you would then need to load the following modules: dealii,
boost, mpi, cmake, python, py-mpi4py, py-matplotlib, py-parsing, py-numpy, and
py-h5py.

Message Passing Interface (MPI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cap should be working with any of the MPI implementations. It has only been
tested with Open MPI, MPICH, and Intel MPI.

Boost
^^^^^
Boost version 1.59.0 or later is required. Boost can be downloaded from `here
<http://www.boost.org/users/download>`_. Make sure to install **all** the
libraries. Do not forget to add the ``using mpi ;`` directive to your
``project-config.jam`` file before building.

deal.II
^^^^^^^
The open source finite element library deal.II is optional. It is only required
to work with energy storage devices of type ``SuperCapacitor``. Version 8.4.0 or
later compiled with C++14/MPI/Boost/p4est/Trilinos support is required. The
development sources can be found `here <https://github.com/dealii/dealii>`_.
Please refer to the deal.II documentation to see how to install `p4est
<https://dealii.org/developer/external-libs/p4est.html>`_ and `Trilinos
<https://dealii.org/developer/external-libs/trilinos.html>`_.


Install Cap from source
-----------------------
Get the source:

.. code::

    $ git clone https://github.com/ORNL-CEES/Cap.git && cd Cap

Create a ``configure_cap.sh`` script such as:

.. code-block:: bash
    :linenos:

    #!/usr/bin/env bash

    EXTRA_ARGS=$@

    cmake \
        -D CMAKE_INSTALL_PREFIX=<your/install/prefix/here> \
        -D BOOST_DIR=<path/to/boost> \
        -D ENABLE_DEAL_II=ON \
        -D DEAL_II_DIR=<path/to/dealii> \
        $EXTRA_ARGS \ 
        ..

Configure, build and install:

.. code::

    $ mkdir build && cd build
    $ vi configure_cap.sh
    $ chmod +x configure_cap.sh
    $ ./configure_cap.sh
    $ make -j<N> && make install


Run the tests:

.. code::

    $ ctest -j<N>


Enable the Python wrappers
--------------------------

To build the Python wrappers Cap must be configured with an extra flag ``-D
ENABLE_PYTHON=ON``. It is recommended to use Python 3.X but Cap has been
successfully built with Python 2.X in the past.

.. code::

    $ ../configure_cap.sh -D ENABLE_PYTHON=ON
    $ make install

Prepend the Python install directory to your ``PYTHONPATH`` environment variable
in order to import the pycap module from your Python interpreter.

.. code::

    $ export PYTHONPATH=<cap/install/prefix>/lib/pythonX.Y/site-packages:${PYTHONPATH}

``X.Y`` stands for the version of Python that was used to build Cap, 
for example 2.7 or 3.5.

Launch Python and try:

.. testcode::

    >>> import pycap
    >>> help(pycap)

Note that a number of Python packages are required to use pycap: numpy,
matplotlib, mpi4py, and h5py.


Build this documentation
------------------------

Run the configuration script with the extra flag:

.. code::

    $ ../configure_cap.sh -D ENABLE_DOCUMENTATION=ON

Open the file ``index.html`` in the directory ``docs/html``.

