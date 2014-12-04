INCLUDE(${PROJECT_SOURCE_DIR}/cmake/macros.cmake)
ADD_CAP_TEST(test_parse_params)
ADD_CAP_TEST(test_constant_current_charge_constant_voltage_discharge)
ADD_CAP_TEST(test_constant_current_cycling)

FILE(COPY ${CAP_DATA_DIR}/meshes/cap_2d_4.ucd DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(RENAME ${CMAKE_CURRENT_BINARY_DIR}/cap_2d_4.ucd ${CMAKE_CURRENT_BINARY_DIR}/mesh_2d.ucd)
