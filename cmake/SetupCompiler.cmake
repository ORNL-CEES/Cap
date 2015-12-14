#### Message Passing Interface (MPI) #########################################
IF(DEFINED MPI_INSTALL_DIR)
    MESSAGE("MPI_INSTALL_DIR=${MPI_INSTALL_DIR}")
    INCLUDE(SetupMPI)
ELSE()
    MESSAGE(FATAL_ERROR "MPI is required. Reconfigure with '-DMPI_INSTALL_DIR=/PATH/TO/MPI'")
ENDIF()

MESSAGE("Setting MPI_CXX_COMPILER as CMAKE_CXX_COMPILER")
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(--std=c++11 COMPILER_SUPPORTS_CXX11)
IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
    MESSAGE(FATAL_ERROR "Compiler has no C++11 support. Please use a different C++ compiler.")
ENDIF()
