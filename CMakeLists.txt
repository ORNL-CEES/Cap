SET(DEAL_II_DIR /Users/qdi/Scratch/install/dealii/dev)
SET( CAP_DATA_DIR /Users/qdi/Projects/cap )

# Set the name of the project and target:
SET(TARGET "cap")

INCLUDE_DIRECTORIES(include)

SET(BOOST_INSTALL_DIR /Users/qdi/Projects/dtk/testing-clang/install/dbg/boost)
# setup boost
INCLUDE_DIRECTORIES(${BOOST_INSTALL_DIR}/include)
SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
FIND_LIBRARY(BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options PATHS ${BOOST_INSTALL_DIR}/lib NO_DEFAULT_PATH)
IF (NOT BOOST_PROGRAM_OPTIONS_LIB)
    MESSAGE(FATAL_ERROR "Boost program_options library not found in ${BOOST_INSTALL_DIR}/lib")
ELSE()
    MESSAGE("Using boost program_options library ${BOOST_PROGRAM_OPTIONS_LIB}")
ENDIF()

# Declare all source files the target consists of:
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib )
ADD_LIBRARY( cache
    include/cache/mp_values.cc
    include/cache/boundary_values.cc
    include/cache/dof_extractor.cc
    include/cache/operator.cc
    include/cache/thermal_operator.cc
    include/cache/electrochemical_operator.cc
    include/cache/utils.cc
    include/cache/super_capacitor.cc
    )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

MESSAGE(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
#DEAL_II_INVOKE_AUTOPILOT()
#TARGET_LINK_LIBRARIES(${TARGET} ${BOOST_PROGRAM_OPTIONS_LIB})

DEAL_II_SETUP_TARGET( cache ) 

ENABLE_TESTING()
ADD_SUBDIRECTORY( test )
