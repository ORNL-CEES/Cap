SET(PyCap_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/property_tree_wrappers.h
   ${CMAKE_CURRENT_SOURCE_DIR}/energy_storage_device_wrappers.h
)
SET(PyCap_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/property_tree_wrappers.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/energy_storage_device_wrappers.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/export_property_tree.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/export_energy_storage_device.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/python_wrappers.cc
)
set(PyCap_HEADERS ${PyCap_HEADERS} PARENT_SCOPE)
set(PyCap_SOURCES ${PyCap_SOURCES} PARENT_SCOPE)
