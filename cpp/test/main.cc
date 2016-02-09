/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#define BOOST_TEST_NO_MAIN
#include <boost/mpi/environment.hpp>
#include <boost/test/unit_test.hpp>

bool
init_function()
{
    return true;
}

int
main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    return boost::unit_test::unit_test_main(&init_function, argc, argv);
}
