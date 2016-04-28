INCLUDE(${CMAKE_SOURCE_DIR}/cmake/UnitTesting.cmake)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/dummy)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/deal.II/dummy)
SET(CPP_TESTS "")
LIST(APPEND
    CPP_TESTS
    test_parse_params
    test_energy_storage_device
    test_resistor_capacitor_circuit
    test_resistor_capacitor_circuit-2
    test_ragone_chart
    )
IF(ENABLE_DEAL_II)
    LIST(APPEND
        CPP_TESTS
        test_geometry
        test_postprocessor
        test_mp_values
        test_equivalent_circuit
        test_exact_transient_solution
        test_supercapacitor
        )
ENDIF()
FOREACH(TEST_NAME ${CPP_TESTS})
    Cap_ADD_BOOST_TEST(${TEST_NAME})
ENDFOREACH()
# Add tests that are run in parallel
IF(ENABLE_DEAL_II)
  Cap_ADD_BOOST_TEST(test_distributed_energy_storage 1 2 4)
ENDIF()
IF(ENABLE_GSL)
    Cap_ADD_BOOST_TEST(test_impedance_spectroscopy      )
    Cap_ADD_BOOST_TEST(test_discrete_fourier_transform  )
ENDIF()
Cap_COPY_INPUT_FILE(input_ragone_chart.info           cpp/test/data)
Cap_COPY_INPUT_FILE(input_impedance_spectroscopy.info cpp/test/data)
Cap_COPY_INPUT_FILE(input_resistor_capacitor.info     cpp/test/data)
Cap_COPY_INPUT_FILE(series_rc.info                    cpp/test/data)
Cap_COPY_INPUT_FILE(parallel_rc.info                  cpp/test/data)
Cap_COPY_INPUT_FILE(super_capacitor.info              cpp/test/data)
Cap_COPY_INPUT_FILE(verification_problems.info        cpp/test/data)
Cap_COPY_INPUT_FILE(read_mesh.info                    cpp/test/data)
Cap_COPY_INPUT_FILE(generate_mesh.info                cpp/test/data)
Cap_COPY_INPUT_FILE(mesh_2d.ucd                       cpp/test/data)
