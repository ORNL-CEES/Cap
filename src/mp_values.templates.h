#include <cap/mp_values.h>
#include <stdexcept>
#include <tuple>

namespace cap {

template <int dim, int spacedim>
MPValues<dim, spacedim>::
MPValues(MPValuesParameters<dim, spacedim> const & params)
{
    std::shared_ptr<boost::property_tree::ptree const> database = params.database;
    std::shared_ptr<Geometry<dim> const> geometry = params.geometry;
    if (geometry) { // vraiment bidon mais necessaire pour eviter d'appeler le constructeur en bloucle
        for (auto const & m : *(geometry->get_materials()))
        {
            std::string const & material_name = m.first;
            std::vector<dealii::types::material_id> const & material_ids = m.second;
//            std::shared_ptr<boost::property_tree::ptree const> material_database =
//                std::make_shared<boost::property_tree::ptree>(database->get_child(material_name));
            std::shared_ptr<MPValues<dim>> material = buildMaterial<dim>(material_name, database);
            for (dealii::types::material_id id : material_ids)
            {
                auto ret = (this->materials).emplace(id, material);
                if (!ret.second)
                    throw std::runtime_error("material id "+std::to_string(id)+" assigned to multiple times");
            }
        }
    }
}



template <int dim, int spacedim>
void
MPValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           std::vector<double> &        values) const

{
    dealii::types::material_id const material_id = cell->material_id();
    auto got = (this->materials).find(material_id);
    if (got == (this->materials).end())
    {
        std::cout<<material_id<<static_cast<int>(material_id)<<"\n";
        for (auto x : this->materials)
            std::cout<<x.first<<"\n";
        throw std::runtime_error("Invalid material id");
    }
    auto material = got->second;
    material->get_values(key, cell, values);
}

template <int dim, int spacedim>
void
MPValues<dim, spacedim>::
get_values(std::string const &                         key,
           active_cell_iterator const &                cell,
           std::vector<dealii::Tensor<1, spacedim> > & values) const

{
    std::ignore = key;
    std::ignore = cell;
    std::ignore = values;
    throw std::runtime_error("Should have been overloaded..."); 
}

//////////////////////// NEW STUFF ////////////////////////////
template <int dim, int spacedim>
PorousElectrodeMPValues<dim, spacedim>::
PorousElectrodeMPValues(MPValuesParameters<dim, spacedim> const & parameters)
    : NewStuffMPValues<dim, spacedim>(parameters)
{
    std::shared_ptr<boost::property_tree::ptree const> database = parameters.database;
    std::string const material_name = database->get<std::string>("ugly_hack");
    std::shared_ptr<boost::property_tree::ptree const> material_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child(material_name));

    std::string const matrix_phase   = material_database->get<std::string>("matrix_phase"  );
    std::shared_ptr<boost::property_tree::ptree const> matrix_phase_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child(matrix_phase));

    // from the matrix_phase database
    double const differential_capacitance       = matrix_phase_database->get<double>("differential_capacitance"      );
    double const exchange_current_density       = matrix_phase_database->get<double>("exchange_current_density"      );
    double const void_volume_fraction           = matrix_phase_database->get<double>("void_volume_fraction"          );
    double const tortuosity_factor              = matrix_phase_database->get<double>("tortuosity_factor"             );
    double const pores_characteristic_dimension = matrix_phase_database->get<double>("pores_characteristic_dimension");
    double const pores_geometry_factor          = matrix_phase_database->get<double>("pores_geometry_factor"         );
    double const mass_density                   = matrix_phase_database->get<double>("mass_density"                  );
    double const electrical_conductivity        = matrix_phase_database->get<double>("electrical_conductivity"       );
    double const heat_capacity                  = matrix_phase_database->get<double>("heat_capacity"                 );
    double const thermal_conductivity           = matrix_phase_database->get<double>("thermal_conductivity"          );
    double const specific_surface_area_per_unit_volume =
        (1.0 + pores_geometry_factor) * void_volume_fraction / pores_characteristic_dimension;
    (this->properties).emplace("specific_capacitance",
        [specific_surface_area_per_unit_volume, differential_capacitance]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), specific_surface_area_per_unit_volume * differential_capacitance);
        }
        );
    (this->properties).emplace("solid_electrical_conductivity",
        [electrical_conductivity]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), electrical_conductivity);
        }
        );

    // from the solution_phase datatabase
    std::string const solution_phase = material_database->get<std::string>("solution_phase");
    std::shared_ptr<boost::property_tree::ptree const> solution_phase_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child(solution_phase));
    double const electrolyte_conductivity       = solution_phase_database->get<double>("conductivity"                  );
    double const electrolyte_mass_density       = solution_phase_database->get<double>("mass_density"                  );
    (this->properties).emplace("liquid_electrical_conductivity",
        [electrolyte_conductivity, void_volume_fraction, tortuosity_factor]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), void_volume_fraction * electrolyte_conductivity / tortuosity_factor);
        }
        );
    // TODO: not sure where to pull this from
    double const anodic_charge_transfer_coefficient   = database->get<double>("anodic_charge_transfer_coefficient"  ,   0.5         );
    double const cathodic_charge_transfer_coefficient = database->get<double>("cathodic_charge_transfer_coefficient",   0.5         );
    double const faraday_constant                     = database->get<double>("faraday_constant"                    ,   9.64853365e4);
    double const gas_constant                         = database->get<double>("gas_constant"                        ,   8.3144621   );
    double const temperature                          = database->get<double>("temperature"                         , 300.0         );
    (this->properties).emplace("faradaic_reaction_coefficient",
        [specific_surface_area_per_unit_volume, exchange_current_density,
        anodic_charge_transfer_coefficient, cathodic_charge_transfer_coefficient,
        faraday_constant, gas_constant, temperature]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(),
                specific_surface_area_per_unit_volume * exchange_current_density *
                (anodic_charge_transfer_coefficient
                    + cathodic_charge_transfer_coefficient) *
                faraday_constant / (gas_constant * temperature)
                );
        }
        );
    (this->properties).emplace("electron_thermal_voltage",
        [faraday_constant, gas_constant, temperature]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), gas_constant * temperature / faraday_constant);
        }
        );
    std::ignore = heat_capacity;
    std::ignore = thermal_conductivity;
    (this->properties).emplace("density",
        [mass_density, void_volume_fraction, electrolyte_mass_density]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(),
                void_volume_fraction*electrolyte_mass_density+(1.0-void_volume_fraction)*mass_density
                );
        }
        );

}



template <int dim, int spacedim>
MetalFoilMPValues<dim, spacedim>::
MetalFoilMPValues(MPValuesParameters<dim, spacedim> const & parameters)
    : NewStuffMPValues<dim, spacedim>(parameters)
{
    std::shared_ptr<boost::property_tree::ptree const> database = parameters.database;
    std::string const material_name = database->get<std::string>("ugly_hack");
    std::shared_ptr<boost::property_tree::ptree const> material_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child(material_name));

    std::string const metal_foil = material_database->get<std::string>("metal_foil"  );
    std::shared_ptr<boost::property_tree::ptree const> metal_foil_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child(metal_foil));

    double const mass_density                   = metal_foil_database->get<double>("mass_density"                  );
    double const electrical_resistivity         = metal_foil_database->get<double>("electrical_resistivity"        );
    double const heat_capacity                  = metal_foil_database->get<double>("heat_capacity"                 );
    double const thermal_conductivity           = metal_foil_database->get<double>("thermal_conductivity"          );

    auto null_property =
        []
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), 0.0);
        };
    (this->properties).emplace("specific_capacitance"          , null_property);
    (this->properties).emplace("faradaic_reaction_coefficient" , null_property);
    (this->properties).emplace("liquid_electrical_conductivity", null_property);
    (this->properties).emplace("solid_electrical_conductivity",
        [electrical_resistivity]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), 1.0 / electrical_resistivity);
        }
        );
    (this->properties).emplace("density",
        [mass_density]
        (active_cell_iterator const &, std::vector<double> & values)
        {
            std::fill(values.begin(), values.end(), mass_density);
        }
        );
    std::ignore = heat_capacity;
    std::ignore = thermal_conductivity;
}


template <int dim, int spacedim>
NewStuffMPValues<dim, spacedim>::
NewStuffMPValues(MPValuesParameters<dim, spacedim> const & parameters)
    : MPValues<dim, spacedim>(parameters)
{}

template <int dim, int spacedim>
void
NewStuffMPValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           std::vector<double> &        values) const
{
    auto got = (this->properties).find(key);
    if (got == (this->properties).end())
        throw std::runtime_error("Invalid material property "+key);
    auto eval = got->second;
    eval(cell, values);
}



//////////////////////// SUPERCAPACITORS ////////////////////////////
template <int dim, int spacedim>
SuperCapacitorMPValues<dim, spacedim>::
SuperCapacitorMPValues(MPValuesParameters<dim, spacedim> const & parameters)
    : MPValues<dim, spacedim>(parameters)
{
    std::shared_ptr<boost::property_tree::ptree const> database = parameters.database;

    if (!database) throw std::runtime_error("supercapacitor material property values database is empty");
//    SuperCapacitorMPValuesParameters<dim, spacedim> const * super_capacitor_parameters = 
//        dynamic_cast<SuperCapacitorMPValuesParameters<dim, spacedim> const *>(&parameters);
    
    this->separator_material_id            = database->get<dealii::types::material_id>("separator_material_id"        );
    this->anode_electrode_material_id      = database->get<dealii::types::material_id>("anode_electrode_material_id"  );
    this->anode_collector_material_id      = database->get<dealii::types::material_id>("anode_collector_material_id"  );
    this->cathode_electrode_material_id    = database->get<dealii::types::material_id>("cathode_electrode_material_id");
    this->cathode_collector_material_id    = database->get<dealii::types::material_id>("cathode_collector_material_id");

    this->separator_thermal_conductivity   = database->get<double>("separator_thermal_conductivity"  );
    this->electrode_thermal_conductivity   = database->get<double>("electrode_thermal_conductivity"  );
    this->collector_thermal_conductivity   = database->get<double>("collector_thermal_conductivity"  );
    this->separator_density                = database->get<double>("separator_density"               );
    this->electrode_density                = database->get<double>("electrode_density"               );
    this->collector_density                = database->get<double>("collector_density"               );
    this->electrolyte_mass_density         = database->get<double>("electrolyte_mass_density"        );
    this->separator_heat_capacity          = database->get<double>("separator_heat_capacity"         );
    this->electrode_heat_capacity          = database->get<double>("electrode_heat_capacity"         );
    this->collector_heat_capacity          = database->get<double>("collector_heat_capacity"         );

    this->differential_capacitance         = database->get<double>("differential_capacitance"        );
    this->electrode_void_volume_fraction   = database->get<double>("electrode_void_volume_fraction"  );
    this->separator_void_volume_fraction   = database->get<double>("separator_void_volume_fraction"  );
    this->electrolyte_conductivity         = database->get<double>("electrolyte_conductivity"        );
    this->solid_phase_conductivity         = database->get<double>("solid_phase_conductivity"        );
    this->collector_electrical_resistivity = database->get<double>("collector_electrical_resistivity");
    this->separator_tortuosity_factor      = database->get<double>("separator_tortuosity_factor"     );
    this->electrode_tortuosity_factor      = database->get<double>("electrode_tortuosity_factor"     );
    this->bruggemans_coefficient           = database->get<double>("bruggemans_coefficient"          );
    this->pores_characteristic_dimension   = database->get<double>("pores_characteristic_dimension"  );
    this->pores_geometry_factor            = database->get<double>("pores_geometry_factor"           );
    this->specific_surface_area_per_unit_volume =
        (1.0 + this->pores_geometry_factor) * this->electrode_void_volume_fraction / this->pores_characteristic_dimension;
                                                                                                                 
    this->anodic_charge_transfer_coefficient   = database->get<double>("anodic_charge_transfer_coefficient"  ,   0.5         );
    this->cathodic_charge_transfer_coefficient = database->get<double>("cathodic_charge_transfer_coefficient",   0.5         );
    this->faraday_constant                     = database->get<double>("faraday_constant"                    ,   9.64853365e4);
    this->gas_constant                         = database->get<double>("gas_constant"                        ,   8.3144621   );
    this->exchange_current_density             = database->get<double>("exchange_current_density"            ,   0.0         );
    this->temperature                          = database->get<double>("temperature"                         , 300.0         );
}

template <int dim, int spacedim>
void 
SuperCapacitorMPValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           std::vector<double> &        values) const
{
    if (key.compare("thermal_conductivity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), this->separator_thermal_conductivity);
        } else if ((cell->material_id() == this->anode_electrode_material_id) 
            || (cell->material_id() == this->cathode_electrode_material_id)) {
            std::fill(values.begin(), values.end(), this->electrode_thermal_conductivity);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), this->collector_thermal_conductivity);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("density_times_heat_capacity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), this->separator_density*this->separator_heat_capacity);
        } else if ((cell->material_id() == this->anode_electrode_material_id)
            || (cell->material_id() == this->cathode_electrode_material_id)) {
            std::fill(values.begin(), values.end(), this->electrode_density*this->electrode_heat_capacity);
        } else if ((cell->material_id() == this->anode_collector_material_id)
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), this->collector_density*this->collector_heat_capacity);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("density") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), this->separator_void_volume_fraction*this->electrolyte_mass_density+(1.0-this->separator_void_volume_fraction)*this->separator_density);
        } else if ((cell->material_id() == this->anode_electrode_material_id)
            || (cell->material_id() == this->cathode_electrode_material_id)) {
            std::fill(values.begin(), values.end(), this->electrode_void_volume_fraction*this->electrolyte_mass_density+(1.0-this->electrode_void_volume_fraction)*this->electrode_density);
        } else if ((cell->material_id() == this->anode_collector_material_id)
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), this->collector_density);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("specific_capacitance") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if ((cell->material_id() == this->anode_electrode_material_id)
            || (cell->material_id() == this->cathode_electrode_material_id)) {
            std::fill(values.begin(), values.end(), this->specific_surface_area_per_unit_volume*this->differential_capacitance);
        } else if ((cell->material_id() == this->anode_collector_material_id)
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("solid_electrical_conductivity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if ((cell->material_id() == this->anode_electrode_material_id) 
            || (cell->material_id() == this->cathode_electrode_material_id)) {
//            std::fill(values.begin(), values.end(), this->solid_phase_conductivity*std::pow((1.0-this->electrode_void_volume_fraction), this->bruggemans_coefficient));
            std::fill(values.begin(), values.end(), this->solid_phase_conductivity);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 1.0/this->collector_electrical_resistivity);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("liquid_electrical_conductivity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
//            std::fill(values.begin(), values.end(), this->electrolyte_conductivity*std::pow(this->separator_void_volume_fraction, this->bruggemans_coefficient));
            std::fill(values.begin(), values.end(), this->electrolyte_conductivity*this->separator_void_volume_fraction/this->separator_tortuosity_factor);
        } else if ((cell->material_id() == this->anode_electrode_material_id) 
            || (cell->material_id() == this->cathode_electrode_material_id)) {
//            std::fill(values.begin(), values.end(), this->electrolyte_conductivity*std::pow(this->electrode_void_volume_fraction, this->bruggemans_coefficient));
            std::fill(values.begin(), values.end(), this->electrolyte_conductivity*this->electrode_void_volume_fraction/this->electrode_tortuosity_factor);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("faradaic_reaction_coefficient") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if ((cell->material_id() == this->anode_electrode_material_id)
            || (cell->material_id() == this->cathode_electrode_material_id)) {
            std::fill(values.begin(), values.end(),
                this->specific_surface_area_per_unit_volume
                    * this->exchange_current_density * (this->anodic_charge_transfer_coefficient + this->cathodic_charge_transfer_coefficient)
                    * this->faraday_constant / (this->gas_constant * this->temperature)
                );
        } else if ((cell->material_id() == this->anode_collector_material_id)
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("electron_thermal_voltage") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if ((cell->material_id() == this->anode_electrode_material_id)
            || (cell->material_id() == this->cathode_electrode_material_id)) {
            std::fill(values.begin(), values.end(), this->gas_constant * this->temperature / this->faraday_constant);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else {
        throw std::runtime_error("Invalid key '"+key+"' material property");
    } // end if key
}

} // end namespace cap
