set(Cap_VERSION_MAJOR 0)
set(Cap_VERSION_MINOR 1)
set(Cap_VERSION_PATCH 0)
set(Cap_VERSION
    ${Cap_VERSION_MAJOR}.${Cap_VERSION_MINOR}.${Cap_VERSION_PATCH})
message("Cap version: ${Cap_VERSION}")

MESSAGE("Build type: ${CMAKE_BUILD_TYPE}")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
    ADD_DEFINITIONS(-DBOOST_DISABLE_ASSERTS)
ELSEIF(CMAKE_BUILD_TYPE MATCHES "Debug")
    # DO NOTHING
ELSE()
    MESSAGE(FATAL_ERROR
        "Possible values for CMAKE_BUILD_TYPE are Debug and Release"
    )
ENDIF()

include(TrackCapRevisionNumber)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/cmake/CapConfigVersion.cmake
    VERSION ${Cap_VERSION}
    COMPATIBILITY AnyNewerVersion
)
configure_file(cmake/CapConfig.cmake.in
    ${CMAKE_BINARY_DIR}/cmake/CapConfig.cmake
    @ONLY
)
install(
    FILES
        ${CMAKE_BINARY_DIR}/cmake/CapConfig.cmake
        ${CMAKE_BINARY_DIR}/cmake/CapConfigVersion.cmake
    DESTINATION lib/cmake/Cap
)
install(FILES ${CMAKE_SOURCE_DIR}/LICENSE.md DESTINATION ${CMAKE_INSTALL_PREFIX})

