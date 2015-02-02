#include <boost/property_tree/ptree.hpp>

void put_default_parameters(boost::property_tree::ptree & params)
{
    params.put("material_properties.separator_heat_capacity",         1.2528e3 );
    params.put("material_properties.electrode_heat_capacity",         0.93e3   );
    params.put("material_properties.collector_heat_capacity",         2.7e3    );
    params.put("material_properties.separator_density",               3.1404e3 );
    params.put("material_properties.electrode_density",               1.34e3   );
    params.put("material_properties.collector_density",               0.89815e3);
    params.put("material_properties.separator_thermal_conductivity",  0.0019e2 );
    params.put("material_properties.electrode_thermal_conductivity",  0.0011e2 );
    params.put("material_properties.collector_thermal_conductivity",  2.37e2   );

    params.put("material_properties.alpha",                           0.0      );
                                                                     
    params.put("material_properties.differential_capacitance",        6.4e-2   );
    params.put("material_properties.separator_void_volume_fraction",  0.6      );
    params.put("material_properties.electrode_void_volume_fraction",  0.67     );
    params.put("material_properties.electrolyte_conductivity",        0.067    );
    params.put("material_properties.solid_phase_conductivity",       52.1      );
    params.put("material_properties.bruggemans_coefficient",          1.5      );
    params.put("material_properties.pores_characteristic_dimension",  1.5e-9   );
    params.put("material_properties.pores_geometry_factor",           2.0      );

    params.put("boundary_values.charge_potential",                    2.2      );
    params.put("boundary_values.discharge_potential",                 1.4      );
    params.put("boundary_values.charge_current_density",            324.65     );
    params.put("boundary_values.discharge_current_density",        -324.65     );
    params.put("boundary_values.initial_potential",                   1.6      );
                                                                  
    params.put("boundary_values.ambient_temperature",       300.0   );
    params.put("boundary_values.heat_transfer_coefficient",   8.0e-2);

    // USUALLY YOU WILL NOT WANT TO MESS WITH THE THESE ONES
    params.put("solid_potential_component",  0);
    params.put("liquid_potential_component", 1);
    params.put("temperature_component",      2);

    params.put("material_properties.separator_material_id",         2);          
    params.put("material_properties.anode_electrode_material_id",   1);
    params.put("material_properties.anode_collector_material_id",   4);
    params.put("material_properties.cathode_electrode_material_id", 3);      
    params.put("material_properties.cathode_collector_material_id", 5);      

    params.put("boundary_values.separator_material_id",         2);          
    params.put("boundary_values.anode_electrode_material_id",   1);
    params.put("boundary_values.anode_collector_material_id",   4);
    params.put("boundary_values.cathode_electrode_material_id", 3);      
    params.put("boundary_values.cathode_collector_material_id", 5);      

    params.put("boundary_values.cathode_boundary_id",           1);
    params.put("boundary_values.anode_boundary_id",             2);
    params.put("boundary_values.upper_boundary_id",             3);
    params.put("boundary_values.lower_boundary_id",             4);
    params.put("boundary_values.other_boundary_id",             5);

    params.put("geometry.electrode_width", 50.0e-6);
    params.put("geometry.separator_width", 25.0e-6);
    params.put("geometry.collector_width",  5.0e-6);
    params.put("geometry.sandwich_height",  5.0e-6);

}

std::shared_ptr<boost::property_tree::ptree> initialize_database()
{
    std::shared_ptr<boost::property_tree::ptree> database(new boost::property_tree::ptree);
    database->put("verbose",                    true);
    database->put("solid_potential_component",     0);
    database->put("liquid_potential_component",    1);
    database->put("temperature_component",         2);
    database->put("electrochemical_block",         0);
    database->put("thermal_block",                 1);
    database->put("separator_material_id",         2);          
    database->put("anode_electrode_material_id",   1);
    database->put("anode_collector_material_id",   4);
    database->put("cathode_electrode_material_id", 3);      
    database->put("cathode_collector_material_id", 5);      
    database->put("cathode_boundary_id",           1);
    database->put("anode_boundary_id",             2);
    database->put("upper_boundary_id",             3);
    database->put("lower_boundary_id",             4);
    database->put("other_boundary_id",             5);
    return database;
}
