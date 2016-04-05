set(Cap_INCLUDE_DIRS ${Cap_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
set(Cap_INCLUDE_DIRS ${Cap_INCLUDE_DIRS} PARENT_SCOPE)

set(Cap_HEADERS
    ${Cap_HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/supercapacitor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/physics.h
    ${CMAKE_CURRENT_SOURCE_DIR}/electrochemical_physics.h
    PARENT_SCOPE
   )

set(Cap_SOURCES
    ${Cap_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/supercapacitor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/physics.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/electrochemical_physics.cc
    PARENT_SCOPE
   )
