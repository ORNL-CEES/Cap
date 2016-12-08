SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/pycap)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/dummy)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/deal.II/dummy)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/wrappers/dummy)


add_subdirectory(wrappers)

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
include_directories(SYSTEM ${MPI4PY_INCLUDE_DIR})
python_add_module(PyCap ${PyCap_SOURCES})

set(PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

target_link_libraries(PyCap Cap)
target_link_libraries(PyCap ${Boost_PYTHON_LIBRARY})
target_link_libraries(PyCap ${PYTHON_LIBRARIES})
set_target_properties(PyCap PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)
install(TARGETS PyCap DESTINATION ${PYTHON_INSTALL_PREFIX}/pycap)

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
    observer_pattern.py
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
