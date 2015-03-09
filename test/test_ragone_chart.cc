#define BOOST_TEST_MODULE TestRagoneChart
#define BOOST_TEST_MAIN
#include <cap/resistor_capacitor.h>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>



namespace cap {

void scan(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const ratio             = database->get<double>("ratio"            );
    double const initial_voltage   = database->get<double>("initial_voltage"  );
    double const final_voltage     = database->get<double>("final_voltage"    );

    double const power_lower_limit = database->get<double>("power_lower_limit");
    double const power_upper_limit = database->get<double>("power_upper_limit");
    int    const steps             = database->get<int   >("steps"            );

    std::shared_ptr<cap::ParallelRC> rc = std::dynamic_pointer_cast<cap::ParallelRC>(dev);
    for (double power = power_lower_limit; power <= power_upper_limit; power *= ratio)
    {
        dev->reset_current(-power/initial_voltage);
        dev->reset_voltage(initial_voltage);
        // TODO: time step control
        double const C = rc->C;
        double const U = rc->U;
        double const R = rc->R_series;
        double const I = rc->I;
        double time_step = (0.5 * C * U * U) / (power + R * I * I) / steps;
        double voltage = initial_voltage;
        double time = 0.0;
        double energy = 0.0;
        
        int step = 1;
        for ( ; voltage >= final_voltage; time+=time_step, ++step)
        {
            dev->evolve_one_time_step_constant_power(time_step, -power);
            voltage = rc->U;
            energy -= rc->U * rc->I * time_step; // TODO: trapeze
        }
        
        os<<boost::format("  %10.7e  %10.7e  %10.7e  %10d \n")
            % power
            % energy
            % time
            % step
            ;
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_ragone_chart )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_ragone_chart", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // scan the system
    std::fstream fout;
    fout.open("ragone_chart_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> ragone_chart_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("ragone_chart"));
    cap::scan(device, ragone_chart_database, fout);
}    
