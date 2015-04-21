#define BOOST_TEST_MODULE NoName
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/format.hpp>
#include <memory>
#include <iostream>
#include <fstream>

namespace cap {

void charge_and_chill(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const charge_time     = database->get<double>("charge_time"    );
    double const relaxation_time = database->get<double>("relaxation_time");
    double const time_step       = database->get<double>("time_step"      );
    double const charge_current  = database->get<double>("charge_current" );
    for (double time = 0.0 ; time < charge_time; time += time_step)
    {
        dev->evolve_one_time_step_constant_current(time_step, charge_current);
        os<<boost::format("%5.1f  ") % time;
        dev->print_data(os);
    }
    for (double time = charge_time ; time < charge_time+relaxation_time; time += time_step)
    {
        dev->evolve_one_time_step_constant_current(time_step, 0.0);
        os<<boost::format("%5.1f  ") % time;
        dev->print_data(os);
    }
}



void check_sanity(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database)

{
    double const initial_voltage = database->get<double>("initial_voltage");
    double       initial_current;
    double       voltage;
    double       current;
    double const tolerance       = database->get<double>("tolerance"      );
    double const time_step       = database->get<double>("time_step"      );

    dev->reset_voltage(initial_voltage);
    dev->get_current(initial_current);
    dev->get_voltage(voltage);
    BOOST_CHECK_CLOSE(voltage, initial_voltage, tolerance);

    dev->evolve_one_time_step_constant_voltage(time_step, initial_voltage);
    dev->get_current(current);
    BOOST_CHECK_CLOSE(current, initial_current, 1.0e-2); // matches percent_tolerance in NoName::reset_voltage()
    dev->get_voltage(voltage);
    BOOST_CHECK_CLOSE(voltage, initial_voltage, tolerance);

    initial_current = database->get<double>("initial_current");
    dev->evolve_one_time_step_constant_current(time_step, initial_current);
    dev->get_current(current);
    BOOST_CHECK_CLOSE(current, initial_current, 1.0e-2); // mostlikely related to grid resolution.  unsure how much accuracy we should expect.
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_no_name )
{
    // read input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_no_name", *input_database);

    // build an energy storage device
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device = 
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters const>(device_database));

    // charge and then relax
    std::fstream fout;
    fout.open("no_name_data", std::fstream::out);
    std::shared_ptr<boost::property_tree::ptree> operating_conditions_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("operating_conditions"));
    cap::charge_and_chill(device, operating_conditions_database, fout);

    // check sanity
    cap::check_sanity(device, input_database);
}
