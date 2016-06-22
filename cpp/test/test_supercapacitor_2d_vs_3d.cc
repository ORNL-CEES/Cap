/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE Supercapacitor_2d_vs_3d

#include "main.cc"

#include <cap/supercapacitor.h>
#include <cap/utils.h>
#include <boost/test/unit_test.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

BOOST_AUTO_TEST_CASE(test_3d_charge)
{
  // First do a charge on a 2D device and then compare the result with a charge
  // using a 3D device.
  // Parse input file
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);
  boost::property_tree::ptree geometry_database;
  boost::property_tree::info_parser::read_info("generate_mesh.info",
                                               geometry_database);
  // Use the same divisions for 2D and 3D
  std::vector<unsigned int> divisions_2(2);
  // Collector
  divisions_2[0] = 1;
  divisions_2[1] = 3;
  geometry_database.put("collector.divisions", cap::to_string(divisions_2));
  // Anode
  divisions_2[0] = 1;
  divisions_2[1] = 2;
  geometry_database.put("anode.divisions", cap::to_string(divisions_2));
  // Separator
  divisions_2[0] = 1;
  divisions_2[1] = 2;
  geometry_database.put("separator.divisions", cap::to_string(divisions_2));
  // Cathode
  divisions_2[0] = 1;
  divisions_2[1] = 2;
  geometry_database.put("cathode.divisions", cap::to_string(divisions_2));
  // Number of refinements
  geometry_database.put("n_refinements", 2);
  device_database.put_child("geometry", geometry_database);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database,
                                      boost::mpi::communicator());
  double const charge_current = 5e-3;
  double const time_step = 1e-2;
  double voltage_2d;
  for (unsigned int i = 0; i < 5; ++i)
    device->evolve_one_time_step_constant_current(time_step, charge_current);
  device->get_voltage(voltage_2d);

  // Now compute the 3D device
  std::vector<unsigned int> divisions_3(3);
  // Collector
  divisions_3[0] = 1;
  divisions_3[1] = 1;
  divisions_3[2] = 3;
  geometry_database.put("collector.divisions", cap::to_string(divisions_3));
  // Anode
  divisions_3[0] = 1;
  divisions_3[1] = 1;
  divisions_3[2] = 2;
  geometry_database.put("anode.divisions", cap::to_string(divisions_3));
  // Separator
  divisions_3[0] = 1;
  divisions_3[1] = 1;
  divisions_3[2] = 2;
  geometry_database.put("separator.divisions", cap::to_string(divisions_3));
  // Cathode
  divisions_3[0] = 1;
  divisions_3[1] = 1;
  divisions_3[2] = 2;
  geometry_database.put("cathode.divisions", cap::to_string(divisions_3));
  device_database.put_child("geometry", geometry_database);
  device_database.put("dim", 3);

  device = cap::EnergyStorageDevice::build(device_database,
                                           boost::mpi::communicator());
  double voltage_3d;
  for (unsigned int i = 0; i < 5; ++i)
    device->evolve_one_time_step_constant_current(time_step, charge_current);
  device->get_voltage(voltage_3d);

  double const percent_tolerance = 1e-2;
  BOOST_CHECK_CLOSE(voltage_2d, voltage_3d, percent_tolerance);
}
