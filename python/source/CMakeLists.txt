SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/pycap)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/source/dummy)

INCLUDE_DIRECTORIES(../include)

SET(CPP_HEADERS
   property_tree_wrappers.h
   energy_storage_device_wrappers.h
   )
SET(CPP_SOURCES
    property_tree_wrappers.cc
    energy_storage_device_wrappers.cc
    python_wrappers.cc
    )


ADD_LIBRARY(pycap MODULE ${CPP_SOURCES})
SET_TARGET_PROPERTIES(pycap PROPERTIES 
    PREFIX "_"
    SUFFIX ".so"
    )

TARGET_LINK_LIBRARIES(pycap ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} cap)
IF(ENABLE_DEAL_II)
    TARGET_LINK_LIBRARIES(pycap ${DEAL_II_TARGET})
ENDIF()
IF(ENABLE_GSL)
    TARGET_LINK_LIBRARIES(pycap ${GSL_LIBRARIES})
ENDIF()
INSTALL(TARGETS pycap DESTINATION ${CMAKE_INSTALL_PREFIX}/python/pycap)

SET(PYTHON_SOURCES
    __init__.py
    data_helpers.py
    time_evolution.py
    end_criterion.py
    stage.py
    charge_discharge.py
    voltammetry.py
    ragone_plot.py
    impedance_spectroscopy.py
    peak_detection.py
    )
FOREACH(PYTHON_SOURCE ${PYTHON_SOURCES})
    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_BINARY_DIR}/python/pycap/${PYTHON_SOURCE}
        DEPENDS ${CMAKE_SOURCE_DIR}/python/source/${PYTHON_SOURCE}
        COMMAND ${CMAKE_COMMAND} 
        ARGS -E copy ${CMAKE_SOURCE_DIR}/python/source/${PYTHON_SOURCE} ${CMAKE_BINARY_DIR}/python/pycap/${PYTHON_SOURCE}
        COMMENT "Copying ${PYTHON_SOURCE}"
        )
    ADD_CUSTOM_TARGET(
        ${PYTHON_SOURCE} ALL
        DEPENDS ${CMAKE_BINARY_DIR}/python/pycap/${PYTHON_SOURCE}
        )
ENDFOREACH()
INSTALL(
    DIRECTORY ${CMAKE_BINARY_DIR}/python/pycap
    DESTINATION ${CMAKE_INSTALL_PREFIX}/python
    FILES_MATCHING PATTERN "*.py"
    )
