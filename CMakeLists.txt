CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

PROJECT(cap)

IF(NOT DEFINED DEAL_II_INSTALL_DIR)
    MESSAGE(FATAL_ERROR "\n"
        "You must pass a flag -DDEAL_II_INSTALL_DIR=/path/to/your/install\n"
        "to cmake.\n"
        )
ENDIF()
IF(NOT DEFINED CAP_DATA_DIR)
    MESSAGE(FATAL_ERROR "\n"
        "You must pass a flag -DCAP_DATA_DIR=/path/to/dir/with/meshes\n"
        "to cmake.\n"
        )
ENDIF()

INCLUDE_DIRECTORIES(include)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
ADD_LIBRARY(cap
    src/mp_values.cc
    src/boundary_values.cc
    src/dof_extractor.cc
    src/operator.cc
    src/thermal_operator.cc
    src/electrochemical_operator.cc
    src/post_processor.cc
    src/utils.cc
    src/super_capacitor.cc
    )

FIND_PACKAGE(deal.II 8.2 REQUIRED PATHS ${DEAL_II_INSTALL_DIR} NO_DEFAULT_PATH)
INCLUDE(${DEAL_II_TARGET_CONFIG})
INCLUDE_DIRECTORIES(SYSTEM ${DEAL_II_INCLUDE_DIRS})

IF(DEFINED BOOST_INSTALL_DIR)
    MESSAGE("BOOST_INSTALL_DIR = ${BOOST_INSTALL_DIR}")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
    FIND_LIBRARY(Boost_LIBRARIES 
        NAMES boost_unit_test_framework
        PATHS ${BOOST_INSTALL_DIR}/lib
        NO_DEFAULT_PATH
        )
    FIND_PATH(Boost_INCLUDE_DIR boost/version.hpp
        PATHS ${BOOST_INSTALL_DIR}/include
        NO_DEFAULT_PATH
        )
    SET(Boost_LIBRARY_DIR ${BOOST_INSTALL_DIR}/lib)
    MESSAGE("Boost_LIBRARIES = ${Boost_LIBRARIES}")
    MESSAGE("Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
    MESSAGE("Boost_LIBRARY_DIR = ${Boost_LIBRARY_DIR}")
    INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIR})
    LINK_DIRECTORIES(${Boost_LIBRARY_DIR})
    SET(MAKE_TEST_WITH_BOOST ON)
ENDIF()

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)
IF(MAKE_DOC)
    ADD_SUBDIRECTORY(doc)
ENDIF()

INSTALL(TARGETS cap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL(DIRECTORY src/ 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/cap
    FILES_MATCHING 
    PATTERN "*.h"
    PATTERN "*.templates.h" EXCLUDE
)

