#define BOOST_TEST_MODULE TestCyclicVoltammetry
#define BOOST_TEST_MAIN
#include <cap/resistor_capacitor.h>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>



namespace cap {

void print_headers(std::ostream & os)
{
    os<<"# cyclic voltammetry\n";
    os<<boost::format( "# %22s  %22s  %22s  \n")
        % "time_t_[second]"
        % "current_I_[ampere]"
        % "voltage_U_[volt]"
        ;
}



void report(double const time, std::shared_ptr<cap::EnergyStorageDevice const> dev, std::ostream & os = std::cout)
{
    double current;
    double voltage;
    dev->get_current(current);
    dev->get_voltage(voltage);
    os<<boost::format("  %22.15e  %22.15e  %22.15e  \n")
        % time
        % current
        % voltage
        ;
}



void scan(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const scan_rate           = database->get<double>("scan_rate"          );
    double const step_size           = database->get<double>("step_size"          );
    double const initial_voltage     = database->get<double>("initial_voltage"    );
    double const voltage_upper_limit = database->get<double>("voltage_upper_limit");
    double const voltage_lower_limit = database->get<double>("voltage_lower_limit");
    double const final_voltage       = database->get<double>("final_voltage"      );
    int    const cycles              = database->get<int   >("cycles"             );

    double const time_step = step_size / scan_rate;
    double time = 0.0;
    print_headers(os);
    dev->reset_voltage(initial_voltage);
    double voltage = initial_voltage;
    for (int n = 0; n < cycles; ++n)
    {
        for ( ; voltage <= voltage_upper_limit; voltage += step_size, time+=time_step)
        {
            dev->evolve_one_time_step_changing_voltage(time_step, voltage);
            report(time, dev, os);
        }
        for ( ; voltage >= voltage_lower_limit; voltage -= step_size, time+=time_step)
        {
            dev->evolve_one_time_step_changing_voltage(time_step, voltage);
            report(time, dev, os);
        }
        for ( ; voltage <= final_voltage; voltage += step_size, time+=time_step)
        {
            dev->evolve_one_time_step_changing_voltage(time_step, voltage);
            report(time, dev, os);
        }
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_cyclic_voltammetry )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_cyclic_voltammetry", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // scan the system
    std::fstream fout;
    fout.open("cyclic_voltammetry_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> cyclic_voltammetry_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("cyclic_voltammetry"));
    cap::scan(device, cyclic_voltammetry_database, fout);
}    
