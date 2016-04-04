INCLUDE(cmake/CapAddTest.cmake)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/dummy)
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
IF(ENABLE_GSL)
    Cap_ADD_BOOST_TEST(test_impedance_spectroscopy      )
    Cap_ADD_BOOST_TEST(test_discrete_fourier_transform  )
ENDIF()
COPY_CAP_INPUT_FILE(input_ragone_chart.info          )
COPY_CAP_INPUT_FILE(input_impedance_spectroscopy.info)
COPY_CAP_INPUT_FILE(input_resistor_capacitor.info    )
COPY_CAP_INPUT_FILE(series_rc.info                   )
COPY_CAP_INPUT_FILE(parallel_rc.info                 )
COPY_CAP_INPUT_FILE(super_capacitor.info             )
COPY_CAP_INPUT_FILE(verification_problems.info       )

FILE(COPY ${CAP_DATA_DIR}/meshes/cap_2d_4.ucd DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cap_2d_4.ucd ${CMAKE_CURRENT_BINARY_DIR}/mesh_2d.ucd)
FILE(COPY ${CAP_DATA_DIR}/meshes/cap_3d_6.ucd DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cap_3d_6.ucd ${CMAKE_CURRENT_BINARY_DIR}/mesh_3d.ucd)
