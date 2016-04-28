INCLUDE(${CMAKE_SOURCE_DIR}/cmake/UnitTesting.cmake)

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

foreach(TEST_NAME ${PYTHON_TESTS})
    Cap_ADD_PYTHON_TEST(${TEST_NAME})
endforeach()

Cap_COPY_INPUT_FILE(series_rc.info       cpp/test/data)
Cap_COPY_INPUT_FILE(parallel_rc.info     cpp/test/data)
Cap_COPY_INPUT_FILE(super_capacitor.info cpp/test/data)
Cap_COPY_INPUT_FILE(read_mesh.info       cpp/test/data)
Cap_COPY_INPUT_FILE(mesh_2d.ucd          cpp/test/data)
