#define BOOST_TEST_MODULE TestLeakageCurrent
#define BOOST_TEST_MAIN
#include <cap/resistor_capacitor.h>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>



namespace cap {

void report(double const t, std::shared_ptr<cap::EnergyStorageDevice const> dev, std::ostream & os = std::cout)
{
    os<<boost::format("%10.5f  ") % t;
    dev->print_data(os);
}



void measure_direct_leakage_current(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const dc_voltage          = database->get<double>("dc_voltage"         );
    double const initial_voltage     = database->get<double>("initial_voltage"    );
    double const initial_time        = database->get<double>("initial_time"       );
    double const final_time          = database->get<double>("final_time"         );
    double const time_step           = database->get<double>("time_step"          );

    dev->reset_voltage(initial_voltage);
    for (double time = initial_time; time < final_time; time += time_step)
    {
        dev->evolve_one_time_step_constant_voltage(time_step, dc_voltage);
        report(time, dev, os);
    }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_leakage_current )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_leakage_current", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // scan the system
    std::fstream fout;
    fout.open("leakage_current_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> cyclic_voltammetry_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("operating_conditions"));
    cap::measure_direct_leakage_current(device, cyclic_voltammetry_database, fout);
}    
