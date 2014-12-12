
EXTERNALPROJECT_ADD(
    dealii
    GIT_REPOSITORY     https://github.com/dealii/dealii.git
    DOWNLOAD_DIR       ${DEAL_II_SOURCE_DIR}
    SOURCE_DIR         ${DEAL_II_SOURCE_DIR}
    INSTALL_DIR        ${DEAL_II_INSTALL_DIR}
    BUILD_IN_SOURCE    0
)
