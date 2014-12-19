IF(NOT DEFINED DEAL_II_INSTALL_DIR)
    SET(DEAL_II_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/dealii")
ENDIF()

EXTERNALPROJECT_ADD(
    DEAL_II
    GIT_REPOSITORY     https://github.com/dealii/dealii.git
    GIT_TAG            master
#    URL                http://www.ces.clemson.edu/dealii/deal.II-8.1.0.tar.gz
    DOWNLOAD_DIR       ${DEAL_II_SOURCE_DIR}
    SOURCE_DIR         ${DEAL_II_SOURCE_DIR}
    INSTALL_DIR        ${DEAL_II_INSTALL_DIR}
    CMAKE_ARGS         -DCMAKE_INSTALL_PREFIX=${DEAL_II_INSTALL_DIR}
                       -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -D CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    BUILD_IN_SOURCE    0
)

FILE(APPEND "${CMAKE_INSTALL_PREFIX}/TPLs.cmake" "SET(DEAL_II_INSTALL_DIR \"${DEAL_II_INSTALL_DIR}\")\n")
