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
