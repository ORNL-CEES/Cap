EXTERNALPROJECT_ADD(
    cap-data
    GIT_REPOSITORY     https://github.com/dalg24/cap-data.git
    DOWNLOAD_DIR       ${CAP_DATA_SOURCE_DIR}
    SOURCE_DIR         ${CAP_DATA_SOURCE_DIR}
    BUILD_IN_SOURCE    0
    CONFIGURE_COMMAND  ""
    BUILD_COMMAND      ""
    INSTALL_COMMAND    ""
)
