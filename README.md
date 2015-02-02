cap
===
DEPENDENCIES
------------
deal.II version 8.2 compiled with C++11 support. The development sources can be found 
[here](https://github.com/dealii/dealii).

Boost Unit Test Framework (UTF) is required to build all the tests.

CONFIGURE, BUILD and TEST
-------------------------

    $ export PREFIX=/path/to/some/working/directory/for/the/project
    $ mkdir -p ${PREFIX}/source

Install deal.II

    $ cd ${PREFIX}/source
    $ wget https://github.com/dealii/dealii/releases/download/v8.2.1/dealii-8.2.1.tar.gz
    $ tar -xf dealii-8.2.1.tar.gz
    $ mkdir -p ${PREFIX}/build/dealii
    $ cd ${PREFIX}/build/dealii
    $ cmake -DCMAKE_INSTALL_PREFIX=${PREFIX}/install/dealii ${PREFIX}/source/dealii-8.2.1
    $ make install    (alternatively $ make -j<N> install)

Configure and build cap

    $ cd ${PREFIX}
    $ git clone https://github.com/dalg24/cap ${PREFIX}/source/cap
    $ git clone https://github.com/dalg24/cap-data ${PREFIX}/source/cap-data
    $ mkdir build/cap
    $ cd build/cap
    $ cmake -DDEAL_II_INSTALL_DIR=${PREFIX}/install/dealii -DCAP_DATA_DIR=${PREFIX}/source/cap-data ${PREFIX}/source/cap
    $ make -j<N>

Run tests

    $ ctest -j<N>

AUTHORS
-------
* [Damien Lebrun-Grandie](https://github.com/dalg24)
