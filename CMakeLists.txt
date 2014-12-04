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

IF(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
ADD_LIBRARY(cap
    include/cap/mp_values.cc
    include/cap/boundary_values.cc
    include/cap/dof_extractor.cc
    include/cap/operator.cc
    include/cap/thermal_operator.cc
    include/cap/electrochemical_operator.cc
    include/cap/post_processor.cc
    include/cap/utils.cc
    include/cap/super_capacitor.cc
    )


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${DEAL_II_INSTALL_DIR}
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
PROJECT(cap)
#DEAL_II_INVOKE_AUTOPILOT()

DEAL_II_SETUP_TARGET(cap) 

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(doc)

#CTEST_START(Nightly)
#SET(CTEST_DROP_METHOD "http")
#SET(CTEST_DROP_SITE "qdi-imac.ornl.gov")
#SET(CTEST_DROP_LOCATION "/~qdi/CDash/submit.php?project=CAP")
#SET(CTEST_DROP_SITE_CDASH TRUE)
#CTEST_SUBMIT()
