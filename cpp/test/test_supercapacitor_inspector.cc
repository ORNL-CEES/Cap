/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE TestSupercapacitorInsepector

#include "main.cc"

#include <cap/energy_storage_device.h>
#include <cap/supercapacitor.h>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <string>

// Check that the inspector can be used to ouput a distributed mesh

namespace cap
{

BOOST_AUTO_TEST_CASE(test_supercapacitor_inspector)
{
  boost::mpi::communicator comm;
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database, comm);

  cap::SuperCapacitorInspector<2> supercap_inspector;
  supercap_inspector.inspect(device.get());

  // Check that the files exist
  if (comm.rank() == 0)
  {
    for (int i = 0; i < comm.size(); ++i)
      BOOST_TEST(boost::filesystem::exists("solution-0000.000" +
                                           std::to_string(i) + ".vtu"));
    BOOST_TEST(boost::filesystem::exists("solution-0000.pvtu"));
  }

  // Remove the files
  if (comm.rank() == 0)
  {
    for (int i = 0; i < comm.size(); ++i)
    {
      std::string filename("solution-0000.000" + std::to_string(i) + ".vtu");
      // Delete file and check return value error code.
      BOOST_TEST(std::remove(filename.c_str()) == 0);
    }
    BOOST_TEST(std::remove("solution-0000.pvtu") == 0);
  }
}
}
