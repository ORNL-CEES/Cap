cap
===
DEPENDENCIES
------------
| Packages               | Required | Version |
| -----------------------|:--------:| -------:|
| deal.II                | Optional | 8.3.0   |
| Boost                  | Required | 1.59.0  |
| GNU Scientific Library | Optional | 1.16    |
| Python                 | Optional | 2.7     |

deal.II version 8.3 compiled with C++11 support. The development sources can be found 
[here](https://github.com/dealii/dealii).

Boost version at least 1.59.0.  Boost Unit Test Framework (UTF) is required for building all the unit tests
and Boost Python is optional for using the Python wrappers for energy storage devices.

CONFIGURE, BUILD, INSTALL
-------------------------

    $ export PREFIX=/path/to/some/working/directory/for/the/project
    $ mkdir -p ${PREFIX}/source

Install boost

    $ cd ${PREFIX}/source
    $ wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.bz2
    $ tar -xf boost_1_59_0.tar.bz2
    $ cd ${PREFIX}/source/boost_1_59_0
    $ ./bootstrap.sh --prefix=${PREFIX}/install/boost
    $ ./b2 install -j<N> variant=release

Install deal.II

    $ cd ${PREFIX}/source
    $ wget https://github.com/dealii/dealii/releases/download/v8.3.0/dealii-8.3.0.tar.gz
    $ tar -xf dealii-8.3.0.tar.gz
    $ mkdir -p ${PREFIX}/build/dealii
    $ cd ${PREFIX}/build/dealii
    $ cmake -DCMAKE_INSTALL_PREFIX=${PREFIX}/install/dealii -DCMAKE_BUILD_TYPE=Release -DDEAL_II_WITH_CXX11=ON -DBOOST_DIR=${PREFIX}/install/boost ${PREFIX}/source/dealii-8.3.0
    $ make -j<N> install

Configure, build, and install cap

    $ cd ${PREFIX}
    $ git clone https://github.com/dalg24/cap ${PREFIX}/source/cap
    $ git clone https://github.com/dalg24/cap-data ${PREFIX}/source/cap-data
    $ mkdir build/cap
    $ cd build/cap
    $ cmake -DCMAKE_INSTALL_PREFIX=${PREFIX}/install/cap -DBOOST_INSTALL_DIR=${PREFIX}/install/boost -DDEAL_II_INSTALL_DIR=${PREFIX}/install/dealii -DCAP_DATA_DIR=${PREFIX}/source/cap-data ${PREFIX}/source/cap
    $ make -j<N> && make install
    
TEST
----
    
Run the unit tests    
    
    $ ctest -j<N>

PYTHON INTERFACE
----------------

Find out where Python is installed

    $ export PYTHON_INSTALL_DIR=`python -c "import sys; print sys.prefix"`

Configure cap to build the python interface and (re)install

    $ cmake -DPYTHON_INSTALL_DIR=${PYTHON_INSTALL_DIR} ${PREFIX}/source/cap


AUTHORS
-------
* [Damien Lebrun-Grandie](https://github.com/dalg24)
