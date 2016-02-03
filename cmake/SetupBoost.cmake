SET(BOOST_ROOT ${BOOST_INSTALL_DIR})
SET(Boost_NO_SYSTEM_PATHS ON)
FIND_PACKAGE(Boost 1.59.0 REQUIRED
    COMPONENTS 
    python 
    unit_test_framework
    serialization
    mpi
    timer
    system
    chrono
    )
MESSAGE("Boost_LIBRARIES=${Boost_LIBRARIES}")
MESSAGE("Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
