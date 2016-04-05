/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/mp_values.h>
#include <stdexcept>
#include <tuple>

namespace cap
{

template <int dim, int spacedim>
MPValues<dim, spacedim>::MPValues(
    MPValuesParameters<dim, spacedim> const &params)
{
  std::shared_ptr<boost::property_tree::ptree const> database = params.database;
  std::shared_ptr<Geometry<dim> const> geometry = params.geometry;
  if (geometry)
  { // vraiment bidon mais necessaire pour eviter d'appeler le constructeur en
    // bloucle
    for (auto const &m : *(geometry->get_materials()))
    {
      std::string const &material_name = m.first;
      std::vector<dealii::types::material_id> const &material_ids = m.second;
      //            std::shared_ptr<boost::property_tree::ptree const>
      //            material_database =
      //                std::make_shared<boost::property_tree::ptree>(database->get_child(material_name));
      std::shared_ptr<MPValues<dim>> material =
          buildMaterial<dim>(material_name, database);
      for (dealii::types::material_id id : material_ids)
      {
        auto ret = (this->materials).emplace(id, material);
        if (!ret.second)
          throw std::runtime_error("material id " + std::to_string(id) +
                                   " assigned to multiple times");
      }
    }
  }
}

template <int dim, int spacedim>
void MPValues<dim, spacedim>::get_values(std::string const &key,
                                         active_cell_iterator const &cell,
                                         std::vector<double> &values) const

{
  dealii::types::material_id const material_id = cell->material_id();
  auto got = (this->materials).find(material_id);
  if (got == (this->materials).end())
  {
    std::cerr << "failed to find\n";
    std::cerr << "id=" << material_id
              << "  int=" << static_cast<int>(material_id)
              << "  to_string=" << std::to_string(material_id) << "\n";
    std::cerr << "in\n";
    for (auto x : this->materials)
      std::cerr << "id=" << x.first << "  int=" << static_cast<int>(x.first)
                << "  to_string=" << std::to_string(x.first) << "\n";
    throw std::runtime_error("Invalid material id " +
                             std::to_string(material_id));
  }
  auto material = got->second;
  material->get_values(key, cell, values);
}

template <int dim, int spacedim>
void MPValues<dim, spacedim>::get_values(
    std::string const &key, active_cell_iterator const &cell,
    std::vector<dealii::Tensor<1, spacedim>> &values) const

{
  std::ignore = key;
  std::ignore = cell;
  std::ignore = values;
  throw std::runtime_error("Should have been overloaded...");
}

//////////////////////// NEW STUFF ////////////////////////////
auto to_meters = [](double const &cm)
{
  return 1e-2 * cm;
};
auto to_kilograms_per_cubic_meter = [](double const &g_per_cm3)
{
  return 1e3 * g_per_cm3;
};
auto to_farads_per_square_meter = [](double const &uF_per_cm2)
{
  return 1e-2 * uF_per_cm2;
};
auto to_amperes_per_square_meter = [](double const &A_per_cm2)
{
  return 1e4 * A_per_cm2;
};
auto to_ohm_meter = [](double const &o_cm)
{
  return 1e-2 * o_cm;
};

template <int dim, int spacedim>
PorousElectrodeMPValues<dim, spacedim>::PorousElectrodeMPValues(
    MPValuesParameters<dim, spacedim> const &parameters)
    : NewStuffMPValues<dim, spacedim>(parameters)
{
  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters.database;
  std::string const material_name = database->get<std::string>("ugly_hack");
  std::shared_ptr<boost::property_tree::ptree const> material_database =
      std::make_shared<boost::property_tree::ptree>(
          database->get_child(material_name));

  std::string const matrix_phase =
      material_database->get<std::string>("matrix_phase");
  std::shared_ptr<boost::property_tree::ptree const> matrix_phase_database =
      std::make_shared<boost::property_tree::ptree>(
          database->get_child(matrix_phase));

  // from the matrix_phase database
  // clang-format off
  double const differential_capacitance       = to_farads_per_square_meter(matrix_phase_database->get<double>("differential_capacitance"));
  double const exchange_current_density       = to_amperes_per_square_meter(matrix_phase_database->get<double>("exchange_current_density"));
  double const void_volume_fraction           = matrix_phase_database->get<double>("void_volume_fraction");
  double const tortuosity_factor              = matrix_phase_database->get<double>("tortuosity_factor");
  double const pores_characteristic_dimension = to_meters(matrix_phase_database->get<double>("pores_characteristic_dimension"));
  double const pores_geometry_factor          = matrix_phase_database->get<double>("pores_geometry_factor");
  double const mass_density                   = to_kilograms_per_cubic_meter(matrix_phase_database->get<double>("mass_density"));
  double const electrical_conductivity        = 1.0/to_ohm_meter(matrix_phase_database->get<double>("electrical_resistivity"));
  double const heat_capacity                  = matrix_phase_database->get<double>("heat_capacity");
  double const thermal_conductivity           = matrix_phase_database->get<double>("thermal_conductivity");
  // clang-format on
  double const specific_surface_area_per_unit_volume =
      (1.0 + pores_geometry_factor) * void_volume_fraction /
      pores_characteristic_dimension;
  (this->properties)
      .emplace("specific_surface_area",
               [specific_surface_area_per_unit_volume](
                   active_cell_iterator const &, std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           specific_surface_area_per_unit_volume);
               });
  (this->properties)
      .emplace("specific_capacitance",
               [specific_surface_area_per_unit_volume,
                differential_capacitance](active_cell_iterator const &,
                                          std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           specific_surface_area_per_unit_volume *
                               differential_capacitance);
               });
  (this->properties)
      .emplace("solid_electrical_conductivity",
               [electrical_conductivity, void_volume_fraction](
                   active_cell_iterator const &, std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           (1.0 - void_volume_fraction) *
                               electrical_conductivity);
               });

  // from the solution_phase datatabase
  std::string const solution_phase =
      material_database->get<std::string>("solution_phase");
  std::shared_ptr<boost::property_tree::ptree const> solution_phase_database =
      std::make_shared<boost::property_tree::ptree>(
          database->get_child(solution_phase));
  double const electrolyte_conductivity =
      1.0 / to_ohm_meter(
                solution_phase_database->get<double>("electrical_resistivity"));
  double const electrolyte_mass_density = to_kilograms_per_cubic_meter(
      solution_phase_database->get<double>("mass_density"));
  (this->properties)
      .emplace("liquid_electrical_conductivity",
               [electrolyte_conductivity, void_volume_fraction,
                tortuosity_factor](active_cell_iterator const &,
                                   std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           void_volume_fraction * electrolyte_conductivity /
                               tortuosity_factor);
               });
  // TODO: not sure where to pull this from
  // clang-format off
  double const anodic_charge_transfer_coefficient   = database->get<double>("anodic_charge_transfer_coefficient", 0.5);
  double const cathodic_charge_transfer_coefficient = database->get<double>("cathodic_charge_transfer_coefficient", 0.5);
  double const faraday_constant                     = database->get<double>("faraday_constant", 9.64853365e4);
  double const gas_constant                         = database->get<double>("gas_constant", 8.3144621);
  double const temperature                          = database->get<double>("temperature", 300.0);
  // clang-format on
  (this->properties)
      .emplace("faradaic_reaction_coefficient",
               [specific_surface_area_per_unit_volume, exchange_current_density,
                anodic_charge_transfer_coefficient,
                cathodic_charge_transfer_coefficient, faraday_constant,
                gas_constant, temperature](active_cell_iterator const &,
                                           std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           specific_surface_area_per_unit_volume *
                               exchange_current_density *
                               (anodic_charge_transfer_coefficient +
                                cathodic_charge_transfer_coefficient) *
                               faraday_constant / (gas_constant * temperature));
               });
  (this->properties)
      .emplace("electron_thermal_voltage",
               [faraday_constant, gas_constant, temperature](
                   active_cell_iterator const &, std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           gas_constant * temperature / faraday_constant);
               });
  std::ignore = heat_capacity;
  std::ignore = thermal_conductivity;
  (this->properties)
      .emplace("density",
               [mass_density, void_volume_fraction, electrolyte_mass_density](
                   active_cell_iterator const &, std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           void_volume_fraction * electrolyte_mass_density +
                               (1.0 - void_volume_fraction) * mass_density);
               });
  (this->properties)
      .emplace("density_of_active_material",
               [mass_density, void_volume_fraction](
                   active_cell_iterator const &, std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           (1.0 - void_volume_fraction) * mass_density);
               });
}

template <int dim, int spacedim>
MetalFoilMPValues<dim, spacedim>::MetalFoilMPValues(
    MPValuesParameters<dim, spacedim> const &parameters)
    : NewStuffMPValues<dim, spacedim>(parameters)
{
  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters.database;
  std::string const material_name = database->get<std::string>("ugly_hack");
  std::shared_ptr<boost::property_tree::ptree const> material_database =
      std::make_shared<boost::property_tree::ptree>(
          database->get_child(material_name));

  std::string const metal_foil =
      material_database->get<std::string>("metal_foil");
  std::shared_ptr<boost::property_tree::ptree const> metal_foil_database =
      std::make_shared<boost::property_tree::ptree>(
          database->get_child(metal_foil));

  // clang-format off
  double const mass_density           = to_kilograms_per_cubic_meter(metal_foil_database->get<double>("mass_density"));
  double const electrical_resistivity = to_ohm_meter(metal_foil_database->get<double>("electrical_resistivity"));
  double const heat_capacity          = metal_foil_database->get<double>("heat_capacity");
  double const thermal_conductivity   = metal_foil_database->get<double>("thermal_conductivity");
  // clang-format on

  auto null_property = [](active_cell_iterator const &,
                          std::vector<double> &values)
  {
    std::fill(values.begin(), values.end(), 0.0);
  };
  (this->properties).emplace("density_of_active_material", null_property);
  (this->properties).emplace("specific_surface_area", null_property);
  (this->properties).emplace("specific_capacitance", null_property);
  (this->properties).emplace("faradaic_reaction_coefficient", null_property);
  (this->properties).emplace("liquid_electrical_conductivity", null_property);
  (this->properties)
      .emplace("solid_electrical_conductivity",
               [electrical_resistivity](active_cell_iterator const &,
                                        std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(),
                           1.0 / electrical_resistivity);
               });
  (this->properties)
      .emplace("density", [mass_density](active_cell_iterator const &,
                                         std::vector<double> &values)
               {
                 std::fill(values.begin(), values.end(), mass_density);
               });
  std::ignore = heat_capacity;
  std::ignore = thermal_conductivity;
}

template <int dim, int spacedim>
NewStuffMPValues<dim, spacedim>::NewStuffMPValues(
    MPValuesParameters<dim, spacedim> const &parameters)
    : MPValues<dim, spacedim>(parameters)
{
}

template <int dim, int spacedim>
void NewStuffMPValues<dim, spacedim>::get_values(
    std::string const &key, active_cell_iterator const &cell,
    std::vector<double> &values) const
{
  auto got = (this->properties).find(key);
  if (got == (this->properties).end())
    throw std::runtime_error("Invalid material property " + key);
  auto eval = got->second;
  eval(cell, values);
}

} // end namespace cap
