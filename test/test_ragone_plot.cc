#include <cap/super_capacitor.h>
#include <cap/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <string>
#include <vector>

// forward declaration (implementation is in "common.cc")
std::shared_ptr<boost::property_tree::ptree> initialize_database();
void put_default_parameters(boost::property_tree::ptree & params);

int main(int argc, char *argv[])
{
    std::shared_ptr<boost::property_tree::ptree> database = initialize_database();
    database->put("verbose", false);
    cap::SuperCapacitorProblem<2> super_capacitor(database);

    // SETTING PROBLEM PARAMETERS
    std::shared_ptr<boost::property_tree::ptree> in(new boost::property_tree::ptree);
    put_default_parameters(*in);
    in->put("test_case",       3  );
    in->put("time_step",       0.1);
    in->put("initial_time",    0.0);
    in->put("final_time",   1800.0);
    in->put("boundary_values.charge_current_density",     324.65);      
    in->put("boundary_values.discharge_current_density", -324.65);      
    in->put("boundary_values.charge_potential",             2.5 );
    in->put("boundary_values.discharge_potential",          1.4 );
    in->put("boundary_values.ambient_temperature",          0.0 );

    // SOLVING THE PROBLEM
    double time_step = 30.0;
    unsigned int const min_steps = 200;
    unsigned int const max_steps = 400;
    double const max_discharge_current_density = 10000.0;
    double const min_discharge_current_density =    10.0;
    std::ostream & os = std::cout;
    std::shared_ptr<boost::property_tree::ptree> out(new boost::property_tree::ptree);
    for (double discharge_current_density = min_discharge_current_density; 
        discharge_current_density <= max_discharge_current_density; 
        discharge_current_density *= 1.25)
    {

        in->put("time_step", time_step);
        in->put("boundary_values.discharge_current_density", -discharge_current_density);
        super_capacitor.run(in, out);
    double max_temperature = out->get<double>("quantities_of_interest.max_temperature");
    double energy_density  = out->get<double>("quantities_of_interest.energy_density" );
    double power_density   = out->get<double>("quantities_of_interest.power_density"  );
    double efficiency      = out->get<double>("quantities_of_interest.efficiency"     );
    double charge_time     = out->get<double>("quantities_of_interest.charge_time"    );
    double discharge_time  = out->get<double>("discharge_time");
    unsigned int steps     = out->get<unsigned int>("steps");
        if (steps < min_steps)
        {
            time_step = out->get<double>("discharge_time") / max_steps;
            in->put("time_step", time_step);
            super_capacitor.run(in, out);

    max_temperature = out->get<double>("quantities_of_interest.max_temperature");
    energy_density  = out->get<double>("quantities_of_interest.energy_density" );
    power_density   = out->get<double>("quantities_of_interest.power_density"  );
    efficiency      = out->get<double>("quantities_of_interest.efficiency"     );
    charge_time     = out->get<double>("quantities_of_interest.charge_time"    );
    discharge_time  = out->get<double>("discharge_time");
    steps           = out->get<unsigned int>("steps");
        }
        os<<boost::format("  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f   %10d  \n")
            % (energy_density/discharge_time*3600.0)
            % energy_density
            % discharge_time
            % max_temperature
            % discharge_current_density
            % time_step
            % steps
            ;
    }

    return 0;
}

#include "common.cc"
