FIND_PACKAGE(Sphinx REQUIRED)

ADD_CUSTOM_TARGET(docs-html ALL
    ${SPHINX_EXECUTABLE}
        -q -b html
        "${CMAKE_CURRENT_SOURCE_DIR}/source"
        "${CMAKE_CURRENT_BINARY_DIR}/html"
    COMMENT "Building HTML documentation with Sphinx"
    DEPENDS pycap
)

INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/docs"
)

ADD_CUSTOM_TARGET(docs-latex ALL
    ${SPHINX_EXECUTABLE}
        -q -b latex
        "${CMAKE_CURRENT_SOURCE_DIR}/source"
        "${CMAKE_CURRENT_BINARY_DIR}/latex"
    COMMENT "Building LATEX documentation with Sphinx"
    DEPENDS pycap
)

INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/latex"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/docs"
)

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        @ONLY
    )
    ADD_CUSTOM_TARGET(doxygen
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        COMMENT "Generating API documentation with Doxygen"
    )
ENDIF()
