ADD_EXECUTABLE( test_constant_current_charge_constant_voltage_discharge.exe test_constant_current_charge_constant_voltage_discharge.cc )
DEAL_II_SETUP_TARGET( test_constant_current_charge_constant_voltage_discharge.exe )
TARGET_LINK_LIBRARIES( test_constant_current_charge_constant_voltage_discharge.exe cache )
ADD_TEST( "constant_current_charge_constant_voltage_discharge" test_constant_current_charge_constant_voltage_discharge.exe )

ADD_EXECUTABLE( test_constant_current_cycling.exe test_constant_current_cycling.cc )
DEAL_II_SETUP_TARGET( test_constant_current_cycling.exe )
TARGET_LINK_LIBRARIES( test_constant_current_cycling.exe cache )
ADD_TEST( "constant_current_cycling" test_constant_current_cycling.exe )

FILE( COPY ${CAP_DATA_DIR}/meshes/cap_2d_4.ucd DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
FILE( RENAME ${CMAKE_CURRENT_BINARY_DIR}/cap_2d_4.ucd ${CMAKE_CURRENT_BINARY_DIR}/mesh_2d.ucd )

ADD_EXECUTABLE( test_parse_params.exe test_parse_params.cc )
DEAL_II_SETUP_TARGET( test_parse_params.exe )
TARGET_LINK_LIBRARIES( test_parse_params.exe cache )
ADD_TEST( "parse_params" test_parse_params.exe )
