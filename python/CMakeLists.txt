SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/pycap)

INCLUDE_DIRECTORIES(include)
ADD_LIBRARY(pycap MODULE
    source/python_wrappers.cc
    source/property_tree_wrappers.cc
    source/property_tree_wrappers.h
    source/energy_storage_device_wrappers.cc
    source/energy_storage_device_wrappers.h
    )

SET_TARGET_PROPERTIES(pycap PROPERTIES 
    PREFIX "_"
    SUFFIX ".so"
    )

TARGET_LINK_LIBRARIES(pycap ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} cap)
IF(ENABLE_DEAL_II)
    TARGET_LINK_LIBRARIES(pycap deal_II)
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


SET(PYTHON_TESTS
    test_property_tree_wrappers
    test_energy_storage_device_wrappers
    test_end_criterion
    test_time_evolution
    test_stage
    test_charge_discharge
    test_voltammetry
    test_ragone_plot
    test_impedance_spectroscopy
)

FUNCTION(ADD_CAP_PYTHON_TEST TEST_NAME)
    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_NAME}.py
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_NAME}.py ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        COMMENT "Copying ${TEST_NAME}.py"
    )
    ADD_CUSTOM_TARGET(
        ${TEST_NAME}.py ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        )
    ADD_TEST(${TEST_NAME} ${PYTHON_EXECUTABLE} ${TEST_NAME}.py)
ENDFUNCTION()

ENABLE_TESTING()
FOREACH(TEST_NAME ${PYTHON_TESTS})
    ADD_CAP_PYTHON_TEST(${TEST_NAME})
ENDFOREACH()

FUNCTION(COPY_CAP_INPUT_FILE INPUT_FILE PATH_TO_FILE)
    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
        DEPENDS ${CMAKE_SOURCE_DIR}/${PATH_TO_FILE}/${INPUT_FILE}
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_SOURCE_DIR}/${PATH_TO_FILE}/${INPUT_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
        COMMENT "Copying ${INPUT_FILE}"
    )
    ADD_CUSTOM_TARGET(
       ${INPUT_FILE}-py ALL
       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${INPUT_FILE}
    )
ENDFUNCTION()

COPY_CAP_INPUT_FILE(device.xml    python/test/data)
COPY_CAP_INPUT_FILE(eis.xml       python/test/data)
COPY_CAP_INPUT_FILE(series_rc.info       test/data)
COPY_CAP_INPUT_FILE(parallel_rc.info     test/data)
COPY_CAP_INPUT_FILE(super_capacitor.info test/data)
