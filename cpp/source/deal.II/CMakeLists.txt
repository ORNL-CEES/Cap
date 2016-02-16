INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/dummy)

SET(Cap_SOURCES
    ${Cap_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/new_supercapacitor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/physics.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/electrochemical_physics.cc
    PARENT_SCOPE
   )
