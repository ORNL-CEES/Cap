#define BOOST_TEST_MODULE TestRagoneChart
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
void report(double const t, std::shared_ptr<cap::SeriesRC const> rc, std::ostream & os = std::cout)
{
    os<<boost::format("  %10.4f  %10.7f  %10.7f  %10.7f  \n")
        % t
        % rc->I
        % rc->U
        % rc->U_C
        ;
}



namespace cap {

void scan(std::shared_ptr<cap::SeriesRC> rc, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const ratio           = database->get<double>("ratio"              );
    double const initial_voltage = database->get<double>("initial_voltage"    );
    double const final_voltage   = database->get<double>("final_voltage"      );

    double const start_power     = database->get<double>("power_limit" );
    double const finish_energy   = database->get<double>("energy_limit");
    int    const steps           = database->get<int   >("steps"       );

    double power = start_power;
    for ( ; ; power *= ratio)
    {
        rc->I = -power/initial_voltage;
        rc->reset(initial_voltage);
        // TODO: time step control
        double const guess = 
        double time_step = (0.5 * rc->C * rc->U * rc->U) / (power + rc->R * rc->I * rc->I) / steps;
        double voltage = initial_voltage;
        double time = 0.0;
        double energy = 0.0;
        
        int step = 1;
        for ( ; voltage >= final_voltage; time+=time_step, ++step)
        {
            rc->evolve_one_time_step_constant_power(time_step, -power);
            voltage = rc->U;
            energy -= rc->U * rc->I * time_step; // TODO: trapeze
//            report(time, rc, os);
        }
        
        os<<boost::format("  %10.7f  %10.7f  %10.7f  %10d \n")
            % power
            % energy
            % time
            % step
            ;
        if (energy <= finish_energy)
            break;
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_ragone_chart )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_ragone_chart", *database);

    // build an energy storage system
    double const C          = database->get<double>("capacitance"        );
    double const R_series   = database->get<double>("series_resistance"  );

    std::shared_ptr<cap::SeriesRC> rc =
        std::make_shared<cap::SeriesRC>(R_series, C);

    // scan the system
    std::fstream fout;
    fout.open("ragone_chart_data", std::fstream::out);

    cap::scan(rc, database, fout);
}    
