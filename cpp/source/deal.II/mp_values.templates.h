/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/mp_values.h>
#include <cap/utils.h>
#include <deal.II/base/function_parser.h>
#include <random>
#include <tuple>
#include <stdexcept>

namespace cap
{

//////////////////////// COMPOSITE MAT /////////////////////////////////////////
template <int dim>
void CompositeMat<dim>::get_values(std::string const &key,
                                   dealii::FEValues<dim> const &fe_values,
                                   std::vector<double> &values) const

{
  auto cell = fe_values.get_cell();
  dealii::types::material_id const material_id = cell->material_id();
  auto got = _materials.find(material_id);
  if (got == _materials.end())
    throw std::runtime_error("Invalid material id " +
                             std::to_string(material_id));
  auto material = got->second;
  material->get_values(key, fe_values, values);
}

//////////////////////// COMPOSITE PRO /////////////////////////////////////////
template <int dim>
void CompositePro<dim>::get_values(std::string const &key,
                                   dealii::FEValues<dim> const &fe_values,
                                   std::vector<double> &values) const
{
  auto got = _properties.find(key);
  if (got == _properties.end())
    throw std::runtime_error("Invalid material property " + key);
  auto property = got->second;
  property->get_values(key, fe_values, values);
}

//////////////////////// UNIFORM CONSTANT //////////////////////////////////////
template <int dim>
UniformConstantMPValues<dim>::UniformConstantMPValues(double const &val)
    : _val(val)
{
}

template <int dim>
void UniformConstantMPValues<dim>::get_values(
    std::string const &key, dealii::FEValues<dim> const &fe_values,
    std::vector<double> &values) const
{
  std::ignore = key;
  std::ignore = fe_values;
  BOOST_ASSERT_MSG(fe_values.get_quadrature().size() == values.size(),
                   "values must be the same size as the quadrature rule"
                   "in MPValues::get_values()");
  std::fill(values.begin(), values.end(), _val);
}

//////////////////////// SUPERCAPACITOR ////////////////////////////////////////
namespace internal
{
// helper function to build parameter distributions
std::function<double(std::default_random_engine &)>
build_parameter(boost::property_tree::ptree const &parameter_database)
{
  auto const distribution_type =
      parameter_database.get<std::string>("distribution_type");
  if (distribution_type.compare("uniform") == 0)
  {
    auto const range =
        to_vector<double>(parameter_database.get<std::string>("range"));
    BOOST_ASSERT_MSG(range.size() == 2,
                     "Invalid range for constructing an uniform distribution");
    std::uniform_real_distribution<double> distribution(range[0], range[1]);
    // ``mutable`` is required to allow the body to modify the distribution
    // because the operator() is non-const.
    return
        [distribution](std::default_random_engine &generator) mutable -> double
    {
      return distribution(generator);
    };
  }
  else if (distribution_type.compare("normal") == 0)
  {
    auto const mean = parameter_database.get<double>("mean");
    auto const standard_deviation =
        parameter_database.get<double>("standard_deviation");
    std::normal_distribution<double> distribution(mean, standard_deviation);
    // ``mutable`` is required to allow the body to modify the distribution
    // because the operator() is non-const.
    return
        [distribution](std::default_random_engine &generator) mutable -> double
    {
      return distribution(generator);
    };
  }
  else if (distribution_type.compare("lognormal") == 0)
  {
    auto const location = parameter_database.get<double>("location");
    auto const scale = parameter_database.get<double>("scale");
    std::lognormal_distribution<double> distribution(location, scale);
    // ``mutable`` is required to allow the body to modify the distribution
    // because the operator() is non-const.
    return
        [distribution](std::default_random_engine &generator) mutable -> double
    {
      return distribution(generator);
    };
  }
  else
    throw std::runtime_error("Invalid parameter distribution type " +
                             distribution_type);
}

// helper function used to construct homogeneous and inhomogeneous materials
template <int dim>
std::unique_ptr<MPValues<dim>>
build_material_properties(std::string const &material_name,
                          MPValuesParameters<dim> const &params)
{
  boost::property_tree::ptree const &database = *params.database;
  boost::property_tree::ptree const &material_database =
      database.get_child(material_name);
  std::string const type = material_database.get<std::string>("type");
  if (type.compare("porous_electrode") == 0)
  {
    return std::make_unique<PorousElectrodeMPValues<dim>>(material_name,
                                                          params);
  }
  else if (type.compare("permeable_membrane") == 0)
  {
    std::string const matrix_phase =
        database.get<std::string>(material_name + "." + "matrix_phase");
    // There is no separate class for a permeable membrane, we just reuse the
    // porous electrode one but it expects parameters that do not make sense for
    // the separator. We make a copy of the original ptree and add the missing
    // parameters that are required for the electrode.
    auto dummy_database =
        std::make_shared<boost::property_tree::ptree>(*params.database);
    dummy_database->put(matrix_phase + "." + "differential_capacitance", 0.0);
    dummy_database->put(matrix_phase + "." + "exchange_current_density", 0.0);
    dummy_database->put(matrix_phase + "." + "electrical_resistivity",
                        std::numeric_limits<double>::max());
    MPValuesParameters<dim> dummy_params = params;
    dummy_params.database = dummy_database;
    return std::make_unique<PorousElectrodeMPValues<dim>>(material_name,
                                                          dummy_params);
  }
  else if (type.compare("current_collector") == 0)
  {
    return std::make_unique<MetalFoilMPValues<dim>>(material_name, params);
  }
  else
  {
    throw std::runtime_error("Invalid material type " + type);
  }
}
} // end namespace internal

//////////////////////// INHOMOGENEOUS SUPERCAPACITOR //////////////////////////
template <int dim>
InhomogeneousSuperCapacitorMPValues<dim>::InhomogeneousSuperCapacitorMPValues(
    MPValuesParameters<dim> const &params)
{
  boost::property_tree::ptree const &database = *params.database;
  // Build a map material_id -> material_name
  std::map<dealii::types::material_id, std::string> material_map;
  for (auto const &m : *params.geometry->get_materials())
  {
    std::string const &material_name = m.first;
    auto const &material_ids = m.second;
    boost::property_tree::ptree const &material_database =
        database.get_child(material_name);
    std::string const type = material_database.get<std::string>("type");
    for (auto const &id : material_ids)
    {
      auto ret = material_map.emplace(id, material_name);
      if (!ret.second)
        throw std::runtime_error("material id " + std::to_string(id) +
                                 " assigned multiple times");
    }
  }
  auto const parameters = database.get<int>("parameters");
  // Build a map material_id -> parameters
  // TODO: For simplicity let's perturb all parameters for now
  // The map is only a map parameter path in the ptree -> object that emulates
  // the std random distribution (i.e. takes a generator object as argument and
  // returns a double)
  std::map<std::string, std::function<double(std::default_random_engine &)>>
      parameter_map;
  for (int p = 0; p < parameters; ++p)
  {
    boost::property_tree::ptree const &parameter_database =
        database.get_child("parameter_" + std::to_string(p));
    auto const parameter_path = parameter_database.get<std::string>("path");
    // Check that the parameter does exist
    if (!database.get_optional<double>(parameter_path))
      throw std::runtime_error("Parameter path " + parameter_path +
                               " does not exist in the material database");
    auto ret = parameter_map.emplace(
        parameter_path, internal::build_parameter(parameter_database));
    if (!ret.second)
      throw std::runtime_error("Parameter " + parameter_path +
                               "  is present multiple times");
  }
  // Make a copy of the material database so that the original database is not
  // changed.
  auto perturbed_database =
      std::make_shared<boost::property_tree::ptree>(*params.database);
  MPValuesParameters<dim> perturbed_params = params;
  perturbed_params.database = perturbed_database;
  // Traverse the triangulation and build material properties with perturbed
  // parameters in each cell.
  // Construct a random number generator and use the MPI rank as a seed.
  std::default_random_engine generator(
      params.geometry->get_mpi_communicator().rank());
  auto const &triangulation = *params.geometry->get_triangulation();

  for (auto cell : triangulation.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      // Perturb the parameters in the copy of the database
      for (auto const &x : parameter_map)
        perturbed_database->put(x.first, x.second(generator));
      _map.emplace(cell->id(),
                   internal::build_material_properties(
                       material_map[cell->material_id()], perturbed_params));
    }
  }
}

template <int dim>
void InhomogeneousSuperCapacitorMPValues<dim>::get_values(
    std::string const &key, dealii::FEValues<dim> const &fe_values,
    std::vector<double> &values) const
{
  auto cell = fe_values.get_cell();
  auto got = _map.find(cell->id());
  if (got == _map.end())
    throw std::runtime_error("Invalid cell property " + cell->id().to_string());
  auto property = got->second;
  property->get_values(key, fe_values, values);
}

//////////////////////// SUPERCAPACITOR ////////////////////////////////////////
template <int dim>
std::unique_ptr<MPValues<dim>>
SuperCapacitorMPValuesFactory<dim>::build(MPValuesParameters<dim> const &params)
{
  if (!params.database)
    throw std::runtime_error("Missing database in SuperCapacitorMPValues");
  if (!params.geometry)
    throw std::runtime_error("Missing geometry in SuperCapacitorMPValues");
  if ((*params.database).get("inhomogeneous", false))
    return std::make_unique<InhomogeneousSuperCapacitorMPValues<dim>>(params);
  else
    return std::make_unique<SuperCapacitorMPValues<dim>>(params);
}

template <int dim>
SuperCapacitorMPValues<dim>::SuperCapacitorMPValues(
    MPValuesParameters<dim> const &params)
{
  // Fill the map material_id -> MPValues
  for (auto const &m : *params.geometry->get_materials())
  {
    std::string const &material_name = m.first;
    auto const &material_ids = m.second;
    // Build the adequate material properties
    std::shared_ptr<MPValues<dim>> properties =
        internal::build_material_properties(material_name, params);
    // Assign it to material ids
    for (auto const &id : material_ids)
    {
      auto ret = (this->_materials).emplace(id, properties);
      if (!ret.second)
        throw std::runtime_error("material id " + std::to_string(id) +
                                 " assigned multiple times");
    }
  }
}

//////////////////////// POROUS ELECTRODE //////////////////////////////////////
// helpers to convert units
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

template <int dim>
PorousElectrodeMPValues<dim>::PorousElectrodeMPValues(
    std::string const &material_name, MPValuesParameters<dim> const &parameters)
{
  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters.database;
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
  double const electrical_resistivity         = to_ohm_meter(matrix_phase_database->get<double>("electrical_resistivity"));
  double const electrical_conductivity        = (electrical_resistivity>1e300) ? 0. : 1.0/electrical_resistivity;
  double const heat_capacity                  = matrix_phase_database->get<double>("heat_capacity");
  double const thermal_conductivity           = matrix_phase_database->get<double>("thermal_conductivity");
  // clang-format on
  double const specific_surface_area_per_unit_volume =
      (1.0 + pores_geometry_factor) * void_volume_fraction /
      pores_characteristic_dimension;
  (this->_properties)
      .emplace("specific_surface_area",
               std::make_shared<UniformConstantMPValues<dim>>(
                   specific_surface_area_per_unit_volume));
  (this->_properties)
      .emplace("specific_capacitance",
               std::make_shared<UniformConstantMPValues<dim>>(
                   specific_surface_area_per_unit_volume *
                   differential_capacitance));
  (this->_properties)
      .emplace("solid_electrical_conductivity",
               std::make_shared<UniformConstantMPValues<dim>>(
                   (1.0 - void_volume_fraction) * electrical_conductivity));

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
  (this->_properties)
      .emplace("liquid_electrical_conductivity",
               std::make_shared<UniformConstantMPValues<dim>>(
                   void_volume_fraction * electrolyte_conductivity /
                   tortuosity_factor));
  // TODO: not sure where to pull this from
  // clang-format off
  double const anodic_charge_transfer_coefficient   = database->get<double>("anodic_charge_transfer_coefficient", 0.5);
  double const cathodic_charge_transfer_coefficient = database->get<double>("cathodic_charge_transfer_coefficient", 0.5);
  double const faraday_constant                     = database->get<double>("faraday_constant", 9.64853365e4);
  double const gas_constant                         = database->get<double>("gas_constant", 8.3144621);
  double const temperature                          = database->get<double>("temperature", 300.0);
  // clang-format on
  (this->_properties)
      .emplace("faradaic_reaction_coefficient",
               std::make_shared<UniformConstantMPValues<dim>>(
                   specific_surface_area_per_unit_volume *
                   exchange_current_density *
                   (anodic_charge_transfer_coefficient +
                    cathodic_charge_transfer_coefficient) *
                   faraday_constant / (gas_constant * temperature)));
  (this->_properties)
      .emplace("electron_thermal_voltage",
               std::make_shared<UniformConstantMPValues<dim>>(
                   gas_constant * temperature / faraday_constant));
  std::ignore = heat_capacity;
  std::ignore = thermal_conductivity;
  (this->_properties)
      .emplace("density", std::make_shared<UniformConstantMPValues<dim>>(
                              void_volume_fraction * electrolyte_mass_density +
                              (1.0 - void_volume_fraction) * mass_density));
  (this->_properties)
      .emplace("density_of_active_material",
               std::make_shared<UniformConstantMPValues<dim>>(
                   (1.0 - void_volume_fraction) * mass_density));

  auto custom = matrix_phase_database->get_child_optional(
      "custom_liquid_electrical_conductivity");
  if (custom)
    (this->_properties)["liquid_electrical_conductivity"] =
        std::make_shared<FunctionSpaceMPValues<dim>>(*custom);
}

//////////////////////// METAL FOIL ////////////////////////////////////////////
template <int dim>
MetalFoilMPValues<dim>::MetalFoilMPValues(
    std::string const &material_name, MPValuesParameters<dim> const &parameters)
{
  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters.database;
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

  auto null_property = std::make_shared<UniformConstantMPValues<dim>>(0.0);

  (this->_properties).emplace("density_of_active_material", null_property);
  (this->_properties).emplace("specific_surface_area", null_property);
  (this->_properties).emplace("specific_capacitance", null_property);
  (this->_properties).emplace("faradaic_reaction_coefficient", null_property);
  (this->_properties).emplace("liquid_electrical_conductivity", null_property);
  (this->_properties)
      .emplace("solid_electrical_conductivity",
               std::make_shared<UniformConstantMPValues<dim>>(
                   1.0 / electrical_resistivity));
  (this->_properties)
      .emplace("density",
               std::make_shared<UniformConstantMPValues<dim>>(mass_density));
  std::ignore = heat_capacity;
  std::ignore = thermal_conductivity;
}

template <int dim>
FunctionSpaceMPValues<dim>::FunctionSpaceMPValues(
    std::shared_ptr<dealii::Function<dim> const> const &function)
    : _function(function)
{
}

template <int dim>
FunctionSpaceMPValues<dim>::FunctionSpaceMPValues(
    boost::property_tree::ptree const &ptree)
{
  auto function = std::make_shared<dealii::FunctionParser<dim>>();
  auto const variables =
      ptree.get<std::string>("variables", dim == 2 ? "x,y" : "x,y,z");
  auto const expression = ptree.get<std::string>("expression");
  auto const constants =
      to_map<double>(ptree.get<std::string>("constants", ""));

  function->initialize(variables, expression, constants);
  _function = function;
}

template <int dim>
void FunctionSpaceMPValues<dim>::get_values(
    std::string const &key, dealii::FEValues<dim> const &fe_values,
    std::vector<double> &values) const
{
  std::ignore = key;
  unsigned int n_q_points = fe_values.n_quadrature_points;
  auto const &points = fe_values.get_quadrature_points();
  BOOST_ASSERT_MSG(points.size() == values.size(),
                   "values must be the same size as the quadrature rule"
                   "in MPValues::get_values()");
  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    values[q] = _function->value(points[q]);
  }
}

} // end namespace cap
