
EXTERNALPROJECT_ADD(
    eigen
#    URL                http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2
    URL                ${EIGEN_ARCHIVE}
    DOWNLOAD_DIR       ${EIGEN_SOURCE_DIR}
    SOURCE_DIR         ${EIGEN_SOURCE_DIR}
    INSTALL_DIR        ${EIGEN_INSTALL_DIR}
    CMAKE_ARGS         -DCMAKE_INSTALL_PREFIX=${EIGEN_INSTALL_DIR}
    BUILD_IN_SOURCE    0
)

