include(${CMAKE_SOURCE_DIR}/cmake/UnitTesting.cmake)

add_definitions(-DBOOST_TEST_DYN_LINK)

include_directories(${CMAKE_SOURCE_DIR}/cpp/source/dummy)
include_directories(${CMAKE_SOURCE_DIR}/cpp/source/deal.II/dummy)

set(CPP_TESTS "")
list(APPEND
    CPP_TESTS
    test_parse_params
    test_energy_storage_device
    test_resistor_capacitor_circuit
    test_resistor_capacitor_circuit-2
    test_ragone_chart
    test_timer
    )
if(ENABLE_DEAL_II)
    list(APPEND
        CPP_TESTS
        test_geometry
        test_postprocessor
        test_mp_values
        test_equivalent_circuit
        test_exact_transient_solution
        test_supercapacitor
        )
endif()
foreach(TEST_NAME ${CPP_TESTS})
    Cap_ADD_BOOST_TEST(${TEST_NAME})
endforeach()

# Add tests that are run in parallel
if(ENABLE_DEAL_II)
  Cap_ADD_BOOST_TEST(test_distributed_energy_storage 1 2 4)
  Cap_ADD_BOOST_TEST(test_supercapacitor_inspector 2)
endif()

Cap_COPY_INPUT_FILE(input_ragone_chart.info           cpp/test/data)
Cap_COPY_INPUT_FILE(series_rc.info                    cpp/test/data)
Cap_COPY_INPUT_FILE(parallel_rc.info                  cpp/test/data)
Cap_COPY_INPUT_FILE(super_capacitor.info              cpp/test/data)
Cap_COPY_INPUT_FILE(verification_problems.info        cpp/test/data)
Cap_COPY_INPUT_FILE(read_mesh.info                    cpp/test/data)
Cap_COPY_INPUT_FILE(generate_mesh.info                cpp/test/data)
Cap_COPY_INPUT_FILE(mesh_2d.ucd                       cpp/test/data)
