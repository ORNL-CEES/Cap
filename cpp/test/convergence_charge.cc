/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE ExactTransientSolution
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/mp_values.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

namespace cap
{

void compute_parameters(
    std::shared_ptr<boost::property_tree::ptree const> input_database,
    std::shared_ptr<boost::property_tree::ptree> output_database)
{
  double const cm2_to_m2 = 0.0001;
  double const cm_to_m = 0.01;
  double const cross_sectional_area =
      cm2_to_m2 * input_database->get<double>("geometry.geometric_area");
  double const electrode_width =
      cm_to_m * input_database->get<double>("geometry.anode_electrode_thickness");
  double const separator_width =
      cm_to_m * input_database->get<double>("geometry.separator_thickness");

  // getting the material parameters values
  std::shared_ptr<boost::property_tree::ptree> material_properties_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database->get_child("material_properties"));
  cap::MPValuesParameters<2> mp_values_params(material_properties_database);
  std::shared_ptr<boost::property_tree::ptree> geometry_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database->get_child("geometry"));
  mp_values_params.geometry =
      std::make_shared<cap::DummyGeometry<2>>(geometry_database);
  std::shared_ptr<cap::MPValues<2>> mp_values =
      std::shared_ptr<cap::MPValues<2>>(new cap::MPValues<2>(mp_values_params));
  // build dummy cell itertor and set its material id
  dealii::Triangulation<2> triangulation;
  dealii::GridGenerator::hyper_cube(triangulation);
  dealii::DoFHandler<2> dof_handler(triangulation);
  dealii::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  // electrode
  cell->set_material_id(input_database->get<dealii::types::material_id>(
      "geometry.anode_electrode_material_id"));
  std::vector<double> electrode_solid_electrical_conductivity_values(1);
  std::vector<double> electrode_liquid_electrical_conductivity_values(1);
  std::vector<double> electrode_specific_capacitance_values(1);
  std::vector<double> electrode_exchange_current_density_values(1);
  std::vector<double> electrode_electron_thermal_voltage_values(1);
  mp_values->get_values("solid_electrical_conductivity", cell,
                        electrode_solid_electrical_conductivity_values);
  mp_values->get_values("liquid_electrical_conductivity", cell,
                        electrode_liquid_electrical_conductivity_values);
  mp_values->get_values("specific_capacitance", cell,
                        electrode_specific_capacitance_values);
  mp_values->get_values("faradaic_reaction_coefficient", cell,
                        electrode_exchange_current_density_values);
  mp_values->get_values("electron_thermal_voltage", cell,
                        electrode_electron_thermal_voltage_values);
  if (electrode_exchange_current_density_values[0] == 0.0)
    throw std::runtime_error("test assumes faradaic processes are present, "
                             "exchange_current_density has to be non zero");
  double const total_current = -1.0; // normalized
  double const dimensionless_exchange_current_density =
      electrode_exchange_current_density_values[0] *
      std::pow(electrode_width, 2) *
      (1.0 / electrode_solid_electrical_conductivity_values[0] +
       1.0 / electrode_liquid_electrical_conductivity_values[0]);
  double const dimensionless_current_density =
      total_current * electrode_width /
      electrode_liquid_electrical_conductivity_values[0] /
      electrode_electron_thermal_voltage_values[0];
  double const ratio_of_solution_phase_to_matrix_phase_conductivities =
      electrode_liquid_electrical_conductivity_values[0] /
      electrode_solid_electrical_conductivity_values[0];

  output_database->put("dimensionless_current_density",
                       dimensionless_current_density);
  output_database->put("dimensionless_exchange_current_density",
                       dimensionless_exchange_current_density);
  output_database->put("ratio_of_solution_phase_to_matrix_phase_conductivities",
                       ratio_of_solution_phase_to_matrix_phase_conductivities);

  output_database->put("position_normalization_factor", electrode_width);
  output_database->put(
      "time_normalization_factor",
      electrode_specific_capacitance_values[0] *
          (1.0 / electrode_solid_electrical_conductivity_values[0] +
           1.0 / electrode_liquid_electrical_conductivity_values[0]) *
          std::pow(electrode_width, 2));

  // separator
  cell->set_material_id(input_database->get<dealii::types::material_id>(
      "geometry.separator_material_id"));
  std::vector<double> separator_liquid_electrical_conductivity_values(1);
  mp_values->get_values("liquid_electrical_conductivity", cell,
                        separator_liquid_electrical_conductivity_values);

  double const potential_drop_across_the_separator =
      -total_current * separator_width /
      separator_liquid_electrical_conductivity_values[0];
  double const voltage_normalization_factor =
      electrode_electron_thermal_voltage_values[0];
  output_database->put("potential_drop_across_the_separator",
                       potential_drop_across_the_separator);
  output_database->put("voltage_normalization_factor",
                       voltage_normalization_factor);
  output_database->put("cross_sectional_area", cross_sectional_area);
}

void verification_problem(
    std::shared_ptr<cap::EnergyStorageDevice> dev,
    std::shared_ptr<boost::property_tree::ptree const> database,
    std::ostream &os = std::cout)
{
  double dimensionless_current_density =
      database->get<double>("dimensionless_current_density");
  double dimensionless_exchange_current_density =
      database->get<double>("dimensionless_exchange_current_density");
  double ratio_of_solution_phase_to_matrix_phase_conductivities =
      database->get<double>(
          "ratio_of_solution_phase_to_matrix_phase_conductivities");

  int const infty =
      database->get<int>("terms_in_truncation_of_infinite_series");
  double const pi = std::acos(-1.0);

  auto compute_dimensionless_overpotential =
      [infty, pi, &ratio_of_solution_phase_to_matrix_phase_conductivities,
       &dimensionless_exchange_current_density, &dimensionless_current_density](
          double const dimensionless_time, double const dimensionless_position)
  {
    std::vector<double> coefficients(infty);
    for (int n = 0; n < infty; ++n)
    {
      coefficients[n] =
          (ratio_of_solution_phase_to_matrix_phase_conductivities *
               std::cos(n * pi) +
           1.0) /
          (dimensionless_exchange_current_density +
           std::pow(n, 2) * std::pow(pi, 2)) *
          std::cos(n * pi * dimensionless_position) *
          std::exp(-(std::pow(n, 2) * std::pow(pi, 2) +
                     dimensionless_exchange_current_density) *
                   dimensionless_time);
    }
    return dimensionless_current_density *
               (1.0 + ratio_of_solution_phase_to_matrix_phase_conductivities) *
               std::exp(-dimensionless_exchange_current_density *
                        dimensionless_time) /
               dimensionless_exchange_current_density -
           dimensionless_current_density *
               (std::cosh(std::sqrt(dimensionless_exchange_current_density) *
                          (1.0 - dimensionless_position)) +
                ratio_of_solution_phase_to_matrix_phase_conductivities *
                    std::cosh(
                        std::sqrt(dimensionless_exchange_current_density) *
                        dimensionless_position)) /
               (std::sqrt(dimensionless_exchange_current_density) *
                std::sinh(std::sqrt(dimensionless_exchange_current_density))) +
           2.0 * dimensionless_current_density *
               std::accumulate(&(coefficients[1]), &(coefficients[infty]), 0.0);
  };

  auto compute_dimensionless_potential_drop_across_the_electrode =
      [&compute_dimensionless_overpotential,
       &ratio_of_solution_phase_to_matrix_phase_conductivities,
       &dimensionless_current_density](double const dimensionless_time)
  {
    return (compute_dimensionless_overpotential(dimensionless_time, 0.0) +
            ratio_of_solution_phase_to_matrix_phase_conductivities *
                compute_dimensionless_overpotential(dimensionless_time, 1.0) -
            dimensionless_current_density *
                ratio_of_solution_phase_to_matrix_phase_conductivities) /
           (1.0 + ratio_of_solution_phase_to_matrix_phase_conductivities);
  };

  auto compute_dimensionless_interfacial_current_density =
      [infty, pi, &ratio_of_solution_phase_to_matrix_phase_conductivities,
       &dimensionless_exchange_current_density, &dimensionless_current_density](
          double const dimensionless_time, double const dimensionless_position)
  {
    std::vector<double> coefficients(infty);
    for (int n = 0; n < infty; ++n)
    {
      coefficients[n] =
          (ratio_of_solution_phase_to_matrix_phase_conductivities *
               std::cos(n * pi) +
           1.0) /
          (dimensionless_exchange_current_density +
           std::pow(n, 2) * std::pow(pi, 2)) /
          (ratio_of_solution_phase_to_matrix_phase_conductivities + 1.0) *
          std::pow(n, 2) * std::pow(pi, 2) *
          std::cos(n * pi * dimensionless_position) *
          std::exp(-(std::pow(n, 2) * std::pow(pi, 2) +
                     dimensionless_exchange_current_density) *
                   dimensionless_time);
    }
    return -2.0 * dimensionless_current_density *
               std::accumulate(&(coefficients[1]), &(coefficients[infty]),
                               0.0) -
           dimensionless_current_density *
               std::sqrt(dimensionless_exchange_current_density) *
               (std::cosh(std::sqrt(dimensionless_exchange_current_density) *
                          (1.0 - dimensionless_position)) +
                ratio_of_solution_phase_to_matrix_phase_conductivities *
                    std::cosh(
                        std::sqrt(dimensionless_exchange_current_density) *
                        dimensionless_position)) /
               ((ratio_of_solution_phase_to_matrix_phase_conductivities + 1.0) *
                std::sinh(std::sqrt(dimensionless_exchange_current_density)));
  };

  // exact vs computed
  double const charge_current = database->get<double>("charge_current");
  double const charge_time    = database->get<double>("charge_time");
  double const time_step      = database->get<double>("time_step");
  double const epsilon = time_step * 1.0e-4;
  double const cross_sectional_area =
      database->get<double>("cross_sectional_area");
  double const time_normalization_factor =
      database->get<double>("time_normalization_factor");
  double const voltage_normalization_factor =
      database->get<double>("voltage_normalization_factor");
  double potential_drop_across_the_separator =
      database->get<double>("potential_drop_across_the_separator");
  potential_drop_across_the_separator *= charge_current / cross_sectional_area;
  dimensionless_current_density *= charge_current / cross_sectional_area;

  std::cout << "delta=" << dimensionless_current_density << "\n";
  std::cout << "nu2  =" << dimensionless_exchange_current_density << "\n";
  std::cout << "beta ="
            << ratio_of_solution_phase_to_matrix_phase_conductivities << "\n";
  std::cout << "time step = " << time_step << std::endl;

  double computed_voltage;
  double exact_voltage;
  double const percent_tolerance = database->get<double>("percent_tolerance");
  for (double time = 0.0; time <= charge_time + epsilon; time += time_step)
  {
    double const dimensionless_time =
        (time + time_step) / time_normalization_factor;
    double const dimensionless_potential_drop_across_the_electrode =
        compute_dimensionless_potential_drop_across_the_electrode(
            dimensionless_time);
    exact_voltage = 2.0 * dimensionless_potential_drop_across_the_electrode *
                        voltage_normalization_factor +
                    potential_drop_across_the_separator;
    dev->evolve_one_time_step_constant_current(time_step, charge_current);
    dev->get_voltage(computed_voltage);
    if ((std::abs(time + time_step - 1e-3) < 1e-7) ||
        (std::abs(time + time_step - 2e-3) < 1e-7) ||
        (std::abs(time + time_step - 3e-3) < 1e-7) ||
        (std::abs(time + time_step - 4e-3) < 1e-7) ||
        (std::abs(time + time_step - 5e-3) < 1e-7) ||
        (std::abs(time + time_step - 6e-3) < 1e-7) ||
        (std::abs(time + time_step - 7e-3) < 1e-7) ||
        (std::abs(time + time_step - 8e-3) < 1e-7) ||
        (std::abs(time + time_step - 9e-3) < 1e-7) ||
        (std::abs(time + time_step - 10e-3) < 1e-7))
      os << boost::format("  %22.15e  %22.15e  %22.15e  \n") %
                (time + time_step) % exact_voltage % computed_voltage;
  }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE(test_exact_transient_solution)
{
  // parse input file
  std::shared_ptr<boost::property_tree::ptree> input_database =
      std::make_shared<boost::property_tree::ptree>();
  boost::property_tree::info_parser::read_info("verification_problems.info",
                                               *input_database);

  // build an energy storage system
  std::shared_ptr<boost::property_tree::ptree> device_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database->get_child("device"));
  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(boost::mpi::communicator(),
                                      *device_database);

  // measure discharge curve
  std::fstream fout;
  fout.open("verification_problem_data", std::fstream::out);

  std::shared_ptr<boost::property_tree::ptree> verification_problem_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database->get_child("verification_problem_subramanian"));

  cap::compute_parameters(device_database, verification_problem_database);

  cap::verification_problem(device, verification_problem_database, fout);

  fout.close();
}
