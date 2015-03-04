#define BOOST_TEST_MODULE TestCyclicVoltammetry
#define BOOST_TEST_MAIN
#include <cap/resistor_capacitor.h>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>



// TODO: electrochemical device should have method to report its current state
// so that we can call something like 
// rc->report(time, os);
void report(double const t, std::shared_ptr<cap::ParallelRC const> rc, std::ostream & os = std::cout)
{
    os<<boost::format("  %10.4f  %10.7f  %10.7f  %10.7f  \n")
        % t
        % rc->I
        % rc->U
        % rc->U_C
        ;
}



namespace cap {

void scan(std::shared_ptr<cap::ParallelRC> rc, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const scan_rate           = database->get<double>("scan_rate"          );
    double const step_size           = database->get<double>("step_size"          );
    double const initial_voltage     = database->get<double>("initial_voltage"    );
    double const upper_voltage_limit = database->get<double>("upper_voltage_limit");
    double const lower_voltage_limit = database->get<double>("lower_voltage_limit");
    double const final_voltage       = database->get<double>("final_voltage"      );
    int    const cycles              = database->get<int   >("cycles"             );

    double const time_step = step_size / scan_rate;
    double time = 0.0;
    rc->reset(initial_voltage);
    for (int n = 0; n < cycles; ++n)
    {
        double voltage = initial_voltage;
        for ( ; voltage <= upper_voltage_limit; voltage += step_size, time+=time_step)
        {
            rc->evolve_one_time_step_constant_voltage(time_step, voltage);
            report(time, rc, os);
        }
        for ( ; voltage >= lower_voltage_limit; voltage -= step_size, time+=time_step)
        {
            rc->evolve_one_time_step_constant_voltage(time_step, voltage);
            report(time, rc, os);
        }
        for ( ; voltage <= final_voltage; voltage += step_size, time+=time_step)
        {
            rc->evolve_one_time_step_constant_voltage(time_step, voltage);
            report(time, rc, os);
        }
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_cyclic_voltammetry )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_cyclic_voltammetry", *database);

    // build an energy storage system
    double const C          = database->get<double>("capacitance"        );
    double const R_parallel = database->get<double>("parallel_resistance");
    double const R_series   = database->get<double>("series_resistance"  );

    std::shared_ptr<cap::ParallelRC> rc =
        std::make_shared<cap::ParallelRC>(R_parallel, C);
    rc->R_series = R_series;

    // scan the system
    std::fstream fout;
    fout.open("cyclic_voltammetry_data", std::fstream::out);

    cap::scan(rc, database, fout);
}    
