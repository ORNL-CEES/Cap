INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/dummy)

SET(CAP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/version.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/energy_storage_device.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/resistor_capacitor.cc
    )
IF(ENABLE_DEAL_II)
    SET(CAP_SOURCES
        ${CAP_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/mp_values.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/boundary_values.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/dof_extractor.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/operator.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/thermal_operator.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/electrochemical_operator.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/post_processor.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/geometry.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/supercapacitor.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/equivalent_circuit.cc
        )
ENDIF()
IF(ENABLE_GSL)
    SET(CAP_SOURCES
        ${CAP_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/electrochemical_impedance_spectroscopy.cc
        )
ENDIF()

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
ADD_LIBRARY(cap ${CAP_SOURCES})
IF(ENABLE_DEAL_II)
    DEAL_II_SETUP_TARGET(cap)
ENDIF()

INSTALL(TARGETS cap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL(DIRECTORY .
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/cap
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.templates.h" EXCLUDE
)

