/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE CheckpointRestart

#include "main.cc"

#include <cap/energy_storage_device.h>
#include <cap/resistor_capacitor.h>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>

boost::property_tree::ptree initialize_database()
{
  boost::property_tree::ptree database;
  database.put("series_resistance", 8.);
  database.put("parallel_resistance", 2.);
  database.put("capacitance", 4.);
  return database;
}

boost::property_tree::ptree initialize_zero_database()
{
  boost::property_tree::ptree database;
  database.put("series_resistance", 0.);
  database.put("parallel_resistance", 0.);
  database.put("capacitance", 0.);
  return database;
}

BOOST_AUTO_TEST_CASE(test_series_rc)
{
  boost::mpi::communicator comm;
  std::string filename = "rc_device.txt";
  cap::SeriesRC rc_device(initialize_database(), comm);

  double const charge_current = 5e-3;
  double const time_step = 1e-2;
  for (unsigned int i = 0; i < 3; ++i)
    rc_device.evolve_one_time_step_constant_current(time_step, charge_current);

  // Save the current state
  rc_device.save(filename);

  // Create a new device
  cap::SeriesRC new_rc_device(initialize_zero_database(), comm);
  // Load the previous device
  new_rc_device.load(filename);

  if (comm.rank() == 0)
  {
    // Check the value
    double const tolerance = 1e-6;
    double current;
    double new_current;
    double voltage;
    double new_voltage;
    rc_device.get_current(current);
    new_rc_device.get_current(new_current);
    rc_device.get_voltage(voltage);
    new_rc_device.get_voltage(new_voltage);
    BOOST_CHECK_CLOSE(current, new_current, tolerance);
    BOOST_CHECK_CLOSE(voltage, new_voltage, tolerance);

    // Delete save file
    std::remove(filename.c_str());
  }
}

BOOST_AUTO_TEST_CASE(test_parallel_rc)
{
  boost::mpi::communicator comm;
  std::string filename = "rc_device.txt";
  cap::ParallelRC rc_device(initialize_database(), comm);

  double const charge_current = 5e-3;
  double const time_step = 1e-2;
  for (unsigned int i = 0; i < 3; ++i)
    rc_device.evolve_one_time_step_constant_current(time_step, charge_current);

  // Save the current state
  rc_device.save(filename);

  // Create a new device
  cap::ParallelRC new_rc_device(initialize_zero_database(), comm);
  // Load the preivous device
  new_rc_device.load(filename);

  if (comm.rank() == 0)
  {
    // Check the value
    double const tolerance = 1e-6;
    double current;
    double new_current;
    double voltage;
    double new_voltage;
    rc_device.get_current(current);
    new_rc_device.get_current(new_current);
    rc_device.get_voltage(voltage);
    new_rc_device.get_voltage(new_voltage);
    BOOST_CHECK_CLOSE(current, new_current, tolerance);
    BOOST_CHECK_CLOSE(voltage, new_voltage, tolerance);

    // Delete save file
    std::remove(filename.c_str());
  }
}

BOOST_AUTO_TEST_CASE(test_supercapacitor)
{
  boost::mpi::communicator comm;
  std::string filename = "device.z";
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);
  boost::property_tree::ptree geometry_database;
  boost::property_tree::info_parser::read_info("generate_mesh.info",
                                               geometry_database);

  geometry_database.put("checkpoint", true);
  geometry_database.put("coarse_mesh_filename", "coarse_mesh.z");
  device_database.put_child("geometry", geometry_database);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database, comm);
  double const charge_current = 5e-3;
  double const time_step = 1e-2;
  for (unsigned int i = 0; i < 3; ++i)
    device->evolve_one_time_step_constant_current(time_step, charge_current);
  // Save the current state
  device->save(filename);

  // Create a new device
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);
  geometry_database.clear();
  geometry_database.put("type", "restart");
  geometry_database.put("coarse_mesh_filename", "coarse_mesh.z");
  device_database.put_child("geometry", geometry_database);
  std::shared_ptr<cap::EnergyStorageDevice> new_device =
      cap::EnergyStorageDevice::build(device_database, comm);
  new_device->load(filename);

  // Check that the values where saved and loaded correctly
  double const tolerance = 1e-6;
  double voltage, new_voltage, current, new_current;
  device->get_current(current);
  device->get_voltage(voltage);
  new_device->get_current(new_current);
  new_device->get_voltage(new_voltage);
  BOOST_CHECK_CLOSE(current, new_current, tolerance);
  BOOST_CHECK_CLOSE(voltage, new_voltage, tolerance);

  // Check that we can keep the computation working
  device->evolve_one_time_step_constant_current(time_step, charge_current);
  new_device->evolve_one_time_step_constant_current(time_step, charge_current);
  device->get_current(current);
  device->get_voltage(voltage);
  new_device->get_current(new_current);
  new_device->get_voltage(new_voltage);
  BOOST_CHECK_CLOSE(current, new_current, tolerance);
  BOOST_CHECK_CLOSE(voltage, new_voltage, tolerance);

  // Delete save file
  if (comm.rank() == 0)
    std::remove(filename.c_str());
}
