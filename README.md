cap
===
DEPENDENCIES
------------
deal.II version 8.2 compiled with C++11 support. The development sources can be found 
[here](https://github.com/dealii/dealii).

boost version > 1.49.0.  Boost Unit Test Framework (UTF) is required only for building all the unit tests.

CONFIGURE, BUILD, INSTALL
-------------------------

    $ export PREFIX=/path/to/some/working/directory/for/the/project
    $ mkdir -p ${PREFIX}/source

Install deal.II

    $ cd ${PREFIX}/source
    $ wget https://github.com/dealii/dealii/releases/download/v8.2.1/dealii-8.2.1.tar.gz
    $ tar -xf dealii-8.2.1.tar.gz
    $ mkdir -p ${PREFIX}/build/dealii
    $ cd ${PREFIX}/build/dealii
    $ cmake -DCMAKE_INSTALL_PREFIX=${PREFIX}/install/dealii -DCMAKE_CXX_FLAGS="-std=c++11" ${PREFIX}/source/dealii-8.2.1
    $ make install    (alternatively $ make -j<N> install)

Configure, build, and install cap

    $ cd ${PREFIX}
    $ git clone https://github.com/dalg24/cap ${PREFIX}/source/cap
    $ git clone https://github.com/dalg24/cap-data ${PREFIX}/source/cap-data
    $ mkdir build/cap
    $ cd build/cap
    $ cmake -DCMAKE_INSTALL_PREFIX=${PREFIX}/install/cap -DCMAKE_CXX_FLAGS="-std=c++11" -DDEAL_II_INSTALL_DIR=${PREFIX}/install/dealii -DCAP_DATA_DIR=${PREFIX}/source/cap-data ${PREFIX}/source/cap
    $ make -j<N>
    $ make install
    
TEST
----

Install boost

    $ wget http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.bz2
    $ tar -xf boost_1_57_0.tar.bz2
    $ cd ${PREFIX}/source/boost_1_57_0
    $ ./bootstrap.sh --prefix=${PREFIX}/install/boost
    $ ./b2 install -j<N> variant=release cxxflags="-std=c++11"
    

Configure deal.II with -DBOOST_DIR=${PREFIX}/install/boost and (re)install

Configure cap with boost and build the tests

    $ cd ${PREFIX}/build/cap
    $ cmake -DBOOST_INSTALL_DIR=${PREFIX}/install/boost ${PREFIX}/source/cap
    $ make -j<N>
    
Run the unit tests    
    
    $ ctest -j<N>


AUTHORS
-------
* [Damien Lebrun-Grandie](https://github.com/dalg24)
