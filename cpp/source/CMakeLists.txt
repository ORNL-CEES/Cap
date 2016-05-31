set(Cap_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/version.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/energy_storage_device.h
    ${CMAKE_CURRENT_SOURCE_DIR}/default_inspector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/resistor_capacitor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.h
)
set(Cap_SOURCES
    ${CMAKE_BINARY_DIR}/cpp/source/version.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/energy_storage_device.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/default_inspector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/resistor_capacitor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.cc
)
if(ENABLE_DEAL_II)
    add_subdirectory(deal.II)
endif()
set(Cap_HEADERS ${Cap_HEADERS} PARENT_SCOPE)
set(Cap_SOURCES ${Cap_SOURCES} PARENT_SCOPE)

set(Cap_INCLUDE_DIRS ${Cap_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
set(Cap_INCLUDE_DIRS ${Cap_INCLUDE_DIRS} PARENT_SCOPE)
