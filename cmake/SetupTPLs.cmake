#### cap-data ################################################################
#INCLUDE(ExternalProject)
#INCLUDE(${PROJECT_SOURCE_DIR}/cmake/cap-data.cmake)

IF(NOT DEFINED CAP_DATA_DIR)
    MESSAGE(FATAL_ERROR "\n"
        "You must pass a flag -DCAP_DATA_DIR=/path/to/dir/with/meshes\n"
        "to cmake.\n"
        )
ENDIF()

#### deal.II #################################################################
IF(DEFINED DEAL_II_INSTALL_DIR)
    FIND_PACKAGE(deal.II 8.3 REQUIRED PATHS ${DEAL_II_INSTALL_DIR} NO_DEFAULT_PATH)
    INCLUDE(${DEAL_II_TARGET_CONFIG})
    INCLUDE_DIRECTORIES(SYSTEM ${DEAL_II_INCLUDE_DIRS})
    SET(ENABLE_DEAL_II ON)
    ADD_DEFINITIONS(-DWITH_DEAL_II)
ENDIF()

#### Tasmanian ###############################################################
IF(DEFINED TASMANIAN_INSTALL_DIR)
    MESSAGE("TASMANIAN_INSTALL_DIR = ${TASMANIAN_INSTALL_DIR}")
    INCLUDE(SetupTasmanian)
    SET(ENABLE_TASMANIAN ON)
ENDIF()

#### GNU Scientific Library ##################################################
IF(DEFINED GSL_INSTALL_DIR)
    MESSAGE("GSL_INSTALL_DIR = ${GSL_INSTALL_DIR}")
    INCLUDE(SetupGSL)
    SET(ENABLE_GSL ON)
ENDIF()

#### Boost ###################################################################
IF(DEFINED BOOST_INSTALL_DIR)
    MESSAGE("BOOST_INSTALL_DIR=${BOOST_INSTALL_DIR}")
    INCLUDE(SetupBoost)
ELSE()
    MESSAGE(FATAL_ERROR "Boost is required. Reconfigure with '-DBOOST_INSTALL_DIR=/PATH/TO/BOOST'")
ENDIF()

#### Python ##################################################################
IF(DEFINED PYTHON_INSTALL_DIR)
    MESSAGE("PYTHON_INSTALL_DIR=${PYTHON_INSTALL_DIR}")
    INCLUDE(SetupPython)
    SET(ENABLE_PYTHON ON)
ENDIF()

