/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/mp_values.h>
#include <cap/energy_storage_device.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>
#include <iostream>

namespace cap
{

// reads database for finite element model and write database for equivalent
// circuit model
void compute_equivalent_circuit(
    boost::property_tree::ptree const &input_database,
    boost::property_tree::ptree &output_database)
{
  // TODO: of course we could clear the database or just overwrite but for
  // now let's just throw an exception if it is not empty
  if (!output_database.empty())
    throw std::runtime_error("output_database was not empty...");

  auto to_meters = [](double const &cm)
  {
    return 0.01 * cm;
  };
  auto to_square_meters = [](double const &cm2)
  {
    return 0.0001 * cm2;
  };

  double const cross_sectional_area =
      to_square_meters(input_database.get<double>("geometry.geometric_area"));
  // clang-format off
  double const electrode_width = to_meters(input_database.get<double>("geometry.anode_electrode_thickness"));
  double const separator_width = to_meters(input_database.get<double>("geometry.separator_thickness"      ));
  double const collector_width = to_meters(input_database.get<double>("geometry.anode_collector_thickness"));
  // clang-format on

  // getting the material parameters values
  std::shared_ptr<boost::property_tree::ptree> material_properties_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database.get_child("material_properties"));
  MPValuesParameters<2> mp_values_params(material_properties_database);
  std::shared_ptr<boost::property_tree::ptree> geometry_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database.get_child("geometry"));
  // build dummy cell iterator and set its material id. Because we use a dummy
  // triangulation, we can use MPI_COMM_WORLD.
  std::shared_ptr<dealii::distributed::Triangulation<2>> triangulation(
      new dealii::distributed::Triangulation<2>(MPI_COMM_WORLD));
  dealii::GridGenerator::hyper_cube(*triangulation);
  dealii::FE_Q<2> fe(1);
  dealii::DoFHandler<2> dof_handler(*triangulation);
  dof_handler.distribute_dofs(fe);
  mp_values_params.geometry = std::make_shared<Geometry<2>>(
      triangulation,
      std::make_shared<std::unordered_map<
          std::string, std::set<dealii::types::material_id>>>(
          std::initializer_list<std::pair<
              std::string const, std::set<dealii::types::material_id>>>{
              {"anode", std::set<dealii::types::material_id>{1}},
              {"separator", std::set<dealii::types::material_id>{2}},
              {"cathode", std::set<dealii::types::material_id>{3}},
              {"collector", std::set<dealii::types::material_id>{4, 5}}}),
      std::make_shared<std::unordered_map<
          std::string, std::set<dealii::types::material_id>>>(
          std::initializer_list<std::pair<
              std::string const, std::set<dealii::types::boundary_id>>>{
              {"anode", std::set<dealii::types::boundary_id>{1}},
              {"cathode", std::set<dealii::types::boundary_id>{2}}}));

  std::shared_ptr<MPValues<2>> mp_values =
      std::make_shared<SuperCapacitorMPValues<2>>(mp_values_params);
  dealii::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  dealii::FEValues<2> fe_values(fe, dealii::QGauss<2>(1),
                                dealii::update_default);
  fe_values.reinit(cell);
  // electrode
  cell->set_material_id(1); // <- matches the anode material_id in the
                            //    initializer list
  std::vector<double> electrode_solid_electrical_conductivity_values(1);
  std::vector<double> electrode_liquid_electrical_conductivity_values(1);
  mp_values->get_values("solid_electrical_conductivity", fe_values,
                        electrode_solid_electrical_conductivity_values);
  mp_values->get_values("liquid_electrical_conductivity", fe_values,
                        electrode_liquid_electrical_conductivity_values);
  double const electrode_resistivity =
      (1.0 / electrode_solid_electrical_conductivity_values[0] +
       1.0 / electrode_liquid_electrical_conductivity_values[0] +
       1.0 / (electrode_solid_electrical_conductivity_values[0] +
              electrode_liquid_electrical_conductivity_values[0])) /
      3.0;
  double const electrode_resistance =
      electrode_resistivity * electrode_width / cross_sectional_area;
  std::vector<double> electrode_specific_capacitance_values(1);
  mp_values->get_values("specific_capacitance", fe_values,
                        electrode_specific_capacitance_values);
  double const electrode_capacitance =
      electrode_specific_capacitance_values[0] * electrode_width *
      cross_sectional_area;
  std::vector<double> electrode_exchange_current_density_values(1);
  mp_values->get_values("faradaic_reaction_coefficient", fe_values,
                        electrode_exchange_current_density_values);
  double const electrode_leakage_resistance =
      1.0 / (electrode_exchange_current_density_values[0] * electrode_width *
             cross_sectional_area);
  std::cout << "ELECTRODE\n";
  std::cout << "    specific_capacitance="
            << electrode_specific_capacitance_values[0] << "\n";
  std::cout << "    solid_electrical_conductivity="
            << electrode_solid_electrical_conductivity_values[0] << "\n";
  std::cout << "    liquid_electrical_conductivity="
            << electrode_liquid_electrical_conductivity_values[0] << "\n";
  std::cout << "    exchange_current_density="
            << electrode_exchange_current_density_values[0] << "\n";
  std::cout << "    width=" << electrode_width << "\n";
  std::cout << "    cross_sectional_area=" << cross_sectional_area << "\n";
  // separator
  cell->set_material_id(2); // <- matches the separator material_id in the
                            //    initializer list
  std::vector<double> separator_liquid_electrical_conductivity_values(1);
  mp_values->get_values("liquid_electrical_conductivity", fe_values,
                        separator_liquid_electrical_conductivity_values);
  double const separator_resistivity =
      1.0 / separator_liquid_electrical_conductivity_values[0];
  double const separator_resistance =
      separator_resistivity * separator_width / cross_sectional_area;
  std::cout << "SEPARATOR\n";
  std::cout << "    liquid_electrical_conductivity="
            << separator_liquid_electrical_conductivity_values[0] << "\n";
  std::cout << "    width=" << separator_width << "\n";
  std::cout << "    cross_sectional_area=" << cross_sectional_area << "\n";
  // collector
  cell->set_material_id(4); // <- matches the collector  material_id in the
                            //    initializer list
  std::vector<double> collector_solid_electrical_conductivity_values(1);
  mp_values->get_values("solid_electrical_conductivity", fe_values,
                        collector_solid_electrical_conductivity_values);
  double const collector_resistivity =
      1.0 / collector_solid_electrical_conductivity_values[0];
  double const collector_resistance =
      collector_resistivity * collector_width / cross_sectional_area;
  std::cout << "COLLECTOR\n";
  std::cout << "    solid_electrical_conductivity="
            << collector_solid_electrical_conductivity_values[0] << "\n";
  std::cout << "    width=" << collector_width << "\n";
  std::cout << "    cross_sectional_area=" << cross_sectional_area << "\n";

  std::cout << "electrode_capacitance=" << electrode_capacitance << "\n";
  std::cout << "electrode_resistance=" << electrode_resistance << "\n";
  std::cout << "electrode_leakage_resistance=" << electrode_leakage_resistance
            << "\n";
  std::cout << "separator_resistance=" << separator_resistance << "\n";
  std::cout << "collector_resistance=" << collector_resistance << "\n";

  // compute the effective resistance and capacitance
  double const sandwich_capacitance = electrode_capacitance / 2.0;
  double const sandwich_resistance = 2.0 * electrode_resistance +
                                     separator_resistance +
                                     2.0 * collector_resistance;
  double const sandwich_leakage_resistance = 2.0 * electrode_leakage_resistance;
  std::cout << "sandwich_capacitance=" << sandwich_capacitance << "\n";
  std::cout << "sandwich_resistance=" << sandwich_resistance << "\n";
  std::cout << "sandwich_leakage_resistance=" << sandwich_leakage_resistance
            << "\n";

  output_database.put("capacitance", sandwich_capacitance);
  output_database.put("series_resistance", sandwich_resistance);
  output_database.put("parallel_resistance", sandwich_leakage_resistance);
  if (std::isfinite(sandwich_leakage_resistance))
    output_database.put("type", "ParallelRC");
  else
    output_database.put("type", "SeriesRC");
}

class EquivalentCircuitBuilder : public EnergyStorageDeviceBuilder
{
public:
  EquivalentCircuitBuilder()
  {
    register_energy_storage_device("EquivalentCircuit", this);
  }
  std::unique_ptr<EnergyStorageDevice>
  build(boost::property_tree::ptree const &ptree,
        boost::mpi::communicator const &comm) override
  {
    boost::property_tree::ptree other;
    compute_equivalent_circuit(ptree, other);
    return EnergyStorageDevice::build(other, comm);
  }
} global_EquivalentCircuitBuilder;

} // end namespace cap
