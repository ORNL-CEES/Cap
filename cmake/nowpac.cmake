
EXTERNALPROJECT_ADD(
    nowpac
    GIT_REPOSITORY     git@bitbucket.org:fmaugust/nowpac.git
    GIT_TAG            NOWPACv1.0
    DOWNLOAD_DIR       ${NOWPAC_SOURCE_DIR}
    SOURCE_DIR         ${NOWPAC_SOURCE_DIR}
    INSTALL_DIR        ${NOWPAC_INSTALL_DIR} 
    CMAKE_ARGS         -D EIGEN_INSTALL_DIR=${EIGEN_INSTALL_DIR}
                       -D NLOPT_INSTALL_DIR=${NLOPT_INSTALL_DIR}
    BUILD_IN_SOURCE    0
)
