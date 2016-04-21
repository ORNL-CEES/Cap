/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE EnergyStorageDevice

#include "main.cc"

#include <cap/energy_storage_device.h>
#include <cap/resistor_capacitor.h>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>
#include <boost/mpi/communicator.hpp>
#include <sstream>

// list of valid inputs to build an EnergyStorageDevice
// These are meant as example
std::list<std::string> const valid_device_input = {
    "series_rc.info", "parallel_rc.info",
#ifdef WITH_DEAL_II
    "super_capacitor.info",
#endif
};

BOOST_AUTO_TEST_CASE(test_energy_storage_device_builders)
{
  boost::mpi::communicator world;
  int const size = world.size();
  int const rank = world.rank();
  std::cout << "hello world from rank " << rank << " of size " << size << "\n";
  for (auto const &filename : valid_device_input)
  {
    boost::property_tree::ptree ptree;
    boost::property_tree::info_parser::read_info(filename, ptree);
    BOOST_CHECK_NO_THROW(cap::EnergyStorageDevice::build(world, ptree));
  }

  // invalid type must throw an exception
  boost::property_tree::ptree ptree;
  ptree.put("type", "InvalidDeviceType");
  BOOST_CHECK_THROW(cap::EnergyStorageDevice::build(world, ptree),
                    std::runtime_error);
}

class ExampleInspector : public cap::EnergyStorageDeviceInspector
{
public:
  // get the device type and set the voltage to 1.4 volt
  void inspect(cap::EnergyStorageDevice *device)
  {
    _type = "Unrecognized";
    // use dynamic_cast to find out the actual type
    if (dynamic_cast<cap::SeriesRC *>(device) != nullptr)
      _type = "SeriesRC";
    else if (dynamic_cast<cap::ParallelRC *>(device) != nullptr)
      _type = "ParallelRC";
    // if we make ExampleInspector friend of the derived
    // class for the EnergyStorageDevice we could virtually
    // do anything
    device->evolve_one_time_step_constant_voltage(1.0, 1.4);
  }
  // the type of the device last inspected
  std::string _type;
};

BOOST_AUTO_TEST_CASE(test_energy_storage_device_inspectors)
{
  std::vector<std::string> valid_device_input(2);
  valid_device_input[0] = "series_rc.info";
  valid_device_input[1] = "parallel_rc.info";
  std::vector<std::string> RC_type(2);
  RC_type[0] = "SeriesRC";
  RC_type[1] = "ParallelRC";
  for (unsigned int i = 0; i < 2; ++i)
  {
    boost::mpi::communicator world;
    boost::property_tree::ptree ptree;
    boost::property_tree::info_parser::read_info(valid_device_input[i], ptree);
    auto device = cap::EnergyStorageDevice::build(world, ptree);
    double voltage;
    device->get_voltage(voltage);
    BOOST_TEST(voltage != 1.4);
    ExampleInspector inspector;
    device->inspect(&inspector);
    BOOST_TEST((inspector._type).compare(RC_type[i]) == 0);
    device->get_voltage(voltage);
    BOOST_TEST(voltage == 1.4);
  }
}

// TODO: won't work for SuperCapacitor
#ifdef WITH_DEAL_II
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(test_serialization, 1)
#endif

BOOST_AUTO_TEST_CASE(test_serialization)
{
  for (auto const &filename : valid_device_input)
  {
    boost::property_tree::ptree ptree;
    boost::property_tree::info_parser::read_info(filename, ptree);
    std::shared_ptr<cap::EnergyStorageDevice> original_device =
        cap::EnergyStorageDevice::build(boost::mpi::communicator(), ptree);

    original_device->evolve_one_time_step_constant_voltage(0.1, 2.1);
    double original_voltage;
    double original_current;
    original_device->get_voltage(original_voltage);
    original_device->get_current(original_current);

    try
    {
      std::stringstream ss;
      // save device
      boost::archive::text_oarchive oa(ss);
      oa.register_type<cap::SeriesRC>();
      oa.register_type<cap::ParallelRC>();
      oa << original_device;
      // print the content of the stream to the screen
      std::cout << ss.str() << "\n";
      BOOST_CHECK(!ss.str().empty());

      // restore device
      boost::archive::text_iarchive ia(ss);
      ia.register_type<cap::SeriesRC>();
      ia.register_type<cap::ParallelRC>();
      std::shared_ptr<cap::EnergyStorageDevice> restored_device;
      ia >> restored_device;
      double restored_voltage;
      double restored_current;
      restored_device->get_voltage(restored_voltage);
      restored_device->get_current(restored_current);
      BOOST_CHECK_EQUAL(original_voltage, restored_voltage);
      BOOST_CHECK_EQUAL(original_current, restored_current);
    }
    catch (boost::archive::archive_exception e)
    {
      BOOST_TEST_MESSAGE("unable to serialize the device");
      BOOST_TEST(false);
    }
  }
}
