set(PYTHON_TESTS
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

function(Cap_ADD_PYTHON_TEST TEST_NAME)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
        COMMENT "Copying ${TEST_NAME}.py"
    )
    add_custom_target(
        ${TEST_NAME}.py ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
    )
    if(ARGN)
        set(NUMBER_OF_PROCESSES_TO_EXECUTE ${ARGN})
    else()
        set(NUMBER_OF_PROCESSES_TO_EXECUTE 1)
    endiF()
    foreach(NPROC ${NUMBER_OF_PROCESSES_TO_EXECUTE})
        if(ENABLE_COVERAGE)
            add_test(
                NAME ${TEST_NAME}_py_${NPROC}
                COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${COVERAGE_EXECUTABLE} run --append ${TEST_NAME}.py
            )
        else()
            ADD_TEST(
                NAME ${TEST_NAME}_py_${NPROC}
                COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${PYTHON_EXECUTABLE} ${TEST_NAME}.py
            )
        endif()
        set_tests_properties(${TEST_NAME}_py_${NPROC} PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}")
    endforeach()
endfunction()

foreach(TEST_NAME ${PYTHON_TESTS})
    Cap_ADD_PYTHON_TEST(${TEST_NAME})
endforeach()


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

COPY_CAP_INPUT_FILE(series_rc.info       cpp/test/data)
COPY_CAP_INPUT_FILE(parallel_rc.info     cpp/test/data)
COPY_CAP_INPUT_FILE(super_capacitor.info cpp/test/data)
COPY_CAP_INPUT_FILE(read_mesh.info       cpp/test/data)
COPY_CAP_INPUT_FILE(mesh_2d.ucd          cpp/test/data)
