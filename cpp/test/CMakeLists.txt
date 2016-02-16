INCLUDE(cmake/CapAddTest.cmake)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/cpp/source/dummy)
SET(CPP_TESTS "")
LIST(APPEND
    CPP_TESTS
    test_parse_params
    test_butler_volmer_kinetics
    test_energy_storage_device
    test_resistor_capacitor_circuit
    test_resistor_capacitor_circuit-2
    test_multi_stage
    test_discharge_curve
    test_charge_curve
    test_lissajous_curve
    test_cyclic_voltammetry
    test_leakage_current
    test_ragone_chart
    test_ragone_chart-2
    )
IF(ENABLE_DEAL_II)
    LIST(APPEND
        CPP_TESTS
        test_geometry
        test_postprocessor
        test_mp_values
        test_dirichlet_boundary_operator
        test_equivalent_circuit-2
        test_supercapacitor
        test_exact_transient_solution
        test_exact_transient_solution-2
        )
ENDIF()
FOREACH(TEST_NAME ${CPP_TESTS})
    Cap_ADD_BOOST_TEST(${TEST_NAME})
ENDFOREACH()
IF(ENABLE_GSL)
    Cap_ADD_BOOST_TEST(test_impedance_spectroscopy      )
    Cap_ADD_BOOST_TEST(test_impedance_spectroscopy-2    )
    Cap_ADD_BOOST_TEST(test_discrete_fourier_transform  )
ENDIF()
IF(ENABLE_TASMANIAN)
    Cap_ADD_BOOST_TEST(test_uq_eis                      )
    COPY_CAP_INPUT_FILE(input_uq_eis                )
ENDIF()
COPY_CAP_INPUT_FILE(input_discharge_curve       )
COPY_CAP_INPUT_FILE(input_charge_curve          )
COPY_CAP_INPUT_FILE(input_lissajous_curve       )
COPY_CAP_INPUT_FILE(input_cyclic_voltammetry    )
COPY_CAP_INPUT_FILE(input_leakage_current       )
COPY_CAP_INPUT_FILE(input_ragone_chart          )
COPY_CAP_INPUT_FILE(input_butler_volmer_kinetics)
COPY_CAP_INPUT_FILE(input_multi_stage           )
COPY_CAP_INPUT_FILE(input_impedance_spectroscopy)
COPY_CAP_INPUT_FILE(input_resistor_capacitor    )
COPY_CAP_INPUT_FILE(input_equivalent_circuit    )
COPY_CAP_INPUT_FILE(input_verification_problem  )
COPY_CAP_INPUT_FILE(series_rc.info              )
COPY_CAP_INPUT_FILE(parallel_rc.info            )
COPY_CAP_INPUT_FILE(super_capacitor.info        )
COPY_CAP_INPUT_FILE(equivalent_circuits.info    )
COPY_CAP_INPUT_FILE(verification_problems.info  )

FILE(COPY ${CAP_DATA_DIR}/meshes/cap_2d_4.ucd DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cap_2d_4.ucd ${CMAKE_CURRENT_BINARY_DIR}/mesh_2d.ucd)
FILE(COPY ${CAP_DATA_DIR}/meshes/cap_3d_6.ucd DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cap_3d_6.ucd ${CMAKE_CURRENT_BINARY_DIR}/mesh_3d.ucd)
