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
DEAL_II_INITIALIZE_CACHED_VARIABLES()
DEAL_II_SETUP_TARGET(cap) 

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

