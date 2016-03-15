/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#define BOOST_TEST_MODULE SuperCapacitor:
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/format.hpp>
#include <memory>
#include <iostream>
#include <fstream>

namespace cap {

void check_sanity(std::shared_ptr<cap::EnergyStorageDevice> dev)
{
    for (auto imposed_current : { 10e-3, 5e-3, 2e-3 })
    {
        dev->evolve_one_time_step_constant_current(2.0, imposed_current);
        double measured_current;
        dev->get_current(measured_current);
        BOOST_TEST(imposed_current== measured_current);
    }

    for (auto imposed_voltage : { 1.4, 1.6, 1.8, 2.0, 2.2 })
    {
        dev->evolve_one_time_step_constant_voltage(2.0, imposed_voltage);
        double measured_voltage;
        dev->get_voltage(measured_voltage);
        BOOST_TEST(imposed_voltage == measured_voltage);
    }

    for (auto imposed_power : { 1e-3, 2e-3, })
    {
        dev->evolve_one_time_step_constant_power(2.0, imposed_power);
        double measured_voltage;
        dev->get_voltage(measured_voltage);
        double measured_current;
        dev->get_current(measured_current);
        BOOST_TEST(imposed_power == measured_voltage * measured_current);
    }

    // TODO constant load is not implemented yet
   // for (auto imposed_load : { 100.0, 33.0, })
   // {
   //     dev->evolve_one_time_step_constant_load(2.0, imposed_load);
   //     double measured_voltage;
   //     dev->get_voltage(measured_voltage);
   //     double measured_current;
   //     dev->get_current(measured_current);
   //     BOOST_TEST(imposed_load == - measured_voltage / measured_current);
   // }
}

} // end namespace cap


double constexpr relative_tolerance = 1.0e-2;

BOOST_AUTO_TEST_CASE( test_supercapacitor, *
    boost::unit_test::tolerance(relative_tolerance))
{
    // build an energy storage device
    boost::property_tree::ptree ptree;
    boost::property_tree::info_parser::read_info(
        "super_capacitor.info", ptree);
    ptree.put("type", "New_SuperCapacitor");
    boost::mpi::communicator world; 
    std::shared_ptr<cap::EnergyStorageDevice> supercap =
        cap::buildEnergyStorageDevice(world, ptree);

    // check sanity
    cap::check_sanity(supercap);
}
