SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/pycap)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/dummy)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/dummy)

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

SET(PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
MESSAGE("PYTHON_INSTALL_PREFIX = ${PYTHON_INSTALL_PREFIX}")

TARGET_LINK_LIBRARIES(pycap cap)
TARGET_LINK_LIBRARIES(pycap ${Boost_MPI_LIBRARY})
TARGET_LINK_LIBRARIES(pycap ${Boost_SERIALIZATION_LIBRARY})
TARGET_LINK_LIBRARIES(pycap ${Boost_PYTHON_LIBRARY})
TARGET_LINK_LIBRARIES(pycap ${PYTHON_LIBRARIES})
IF(ENABLE_DEAL_II)
    DEAL_II_SETUP_TARGET(pycap)
ENDIF()
IF(ENABLE_GSL)
    TARGET_LINK_LIBRARIES(pycap ${GSL_LIBRARIES})
ENDIF()
INSTALL(TARGETS pycap DESTINATION ${PYTHON_INSTALL_PREFIX}/pycap)

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
    DESTINATION ${PYTHON_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.py"
    )
