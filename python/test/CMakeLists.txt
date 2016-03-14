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
    test_observer_pattern
)

FUNCTION(ADD_CAP_PYTHON_TEST TEST_NAME)
    ADD_CUSTOM_COMMAND(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        COMMENT "Copying ${TEST_NAME}.py"
    )
    ADD_CUSTOM_TARGET(
        ${TEST_NAME}.py ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
    )
    if(ENABLE_COVERAGE)
        add_test(
            NAME ${TEST_NAME}_py
            COMMAND ${COVERAGE_EXECUTABLE} run --append ${TEST_NAME}.py
        )
    else()
        ADD_TEST(
            NAME ${TEST_NAME}_py
            COMMAND ${PYTHON_EXECUTABLE} ${TEST_NAME}.py
        )
    endif()
ENDFUNCTION()

FOREACH(TEST_NAME ${PYTHON_TESTS})
    ADD_CAP_PYTHON_TEST(${TEST_NAME})
    SET_TESTS_PROPERTIES(${TEST_NAME}_py PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}")
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
COPY_CAP_INPUT_FILE(series_rc.info       cpp/test/data)
COPY_CAP_INPUT_FILE(parallel_rc.info     cpp/test/data)
COPY_CAP_INPUT_FILE(super_capacitor.info cpp/test/data)

