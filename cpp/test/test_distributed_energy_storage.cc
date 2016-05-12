/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE DistributedEnergyStorage

#include "main.cc"

#include <cap/energy_storage_device.h>
#include <cap/mp_values.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

namespace cap
{

void distributed_problem(std::shared_ptr<cap::EnergyStorageDevice> dev)
{
  double const charge_current = 5e-3;
  // This is the values computed using one processor
  double const exact_voltage = 0.24307431815;
  double const time_step = 1e-2;
  double const percent_tolerance = 1e-2;
  double computed_voltage;
  double computed_current;
  for (unsigned int i = 0; i < 3; ++i)
    dev->evolve_one_time_step_constant_current(time_step, charge_current);
  dev->get_current(computed_current);
  dev->get_voltage(computed_voltage);

  BOOST_CHECK_CLOSE(computed_voltage, exact_voltage, percent_tolerance);
  BOOST_CHECK_CLOSE(computed_current, charge_current, percent_tolerance);
}
}

BOOST_AUTO_TEST_CASE(test_distributed_energy_storage)
{
  // Parse input file
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);
  boost::property_tree::ptree geometry_database;
  boost::property_tree::info_parser::read_info("generate_mesh.info",
                                               geometry_database);
  device_database.put_child("geometry", geometry_database);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database,
                                      boost::mpi::communicator());

  cap::distributed_problem(device);
}
