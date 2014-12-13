
EXTERNALPROJECT_ADD(
    dealii
#    GIT_REPOSITORY     https://github.com/dealii/dealii.git
    URL                http://www.ces.clemson.edu/dealii/deal.II-8.1.0.tar.gz
    DOWNLOAD_DIR       ${DEAL_II_SOURCE_DIR}
    SOURCE_DIR         ${DEAL_II_SOURCE_DIR}
    INSTALL_DIR        ${DEAL_II_INSTALL_DIR}
    CMAKE_ARGS         -DCMAKE_INSTALL_PREFIX=${DEAL_II_INSTALL_DIR}
    BUILD_IN_SOURCE    0
)
