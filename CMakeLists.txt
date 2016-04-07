cmake_minimum_required(VERSION 3.1)

project(Cap LANGUAGES CXX VERSION ${Cap_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(SetupTPLs)
include(SetupCap)
if(ENABLE_COVERAGE)
    include(CodeCoverage)
endif()
if(ENABLE_FORMAT)
    include(CodeFormat)
endif()

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

enable_testing()

add_subdirectory(cpp)

if(ENABLE_PYTHON)
    add_subdirectory(python)
endif()

if(ENABLE_DOCUMENTATION)
    add_subdirectory(docs)
endif()

#### OUTDATED ####################################
IF(MAKE_DOC)
    ADD_SUBDIRECTORY(doc-old)
ENDIF()
##################################################
