#include <cap/super_capacitor.h>
#include <cap/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <string>
#include <vector>

// forward declaration (implementation is in "common.cc")
std::shared_ptr<boost::property_tree::ptree> initialize_database();
void put_default_parameters(boost::property_tree::ptree & params);

int main(int argc, char *argv[])
{
    cap::SuperCapacitorProblem<2> super_capacitor(initialize_database());

    // SETTING PROBLEM PARAMETERS
    std::shared_ptr<boost::property_tree::ptree> in(new boost::property_tree::ptree);
    put_default_parameters(*in);
    in->put("test_case",       2  );
    in->put("time_step",       1.0);
    in->put("initial_time",    0.0);
    in->put("final_time",    300.0);
    in->put("max_cycles",    100  );
    in->put("boundary_values.charge_current_density",     324.65);      
    in->put("boundary_values.discharge_current_density", -324.65);      

    // SOLVING THE PROBLEM
    std::shared_ptr<boost::property_tree::ptree> out(new boost::property_tree::ptree);
    super_capacitor.run(in, out);

    // POSTPROCESSING QUANTITIES OF INTEREST
    double max_temperature = out->get<double>("quantities_of_interest.max_temperature");
    double energy_density  = out->get<double>("quantities_of_interest.energy_density" );
    double power_density   = out->get<double>("quantities_of_interest.power_density"  );
    double efficiency      = out->get<double>("quantities_of_interest.efficiency"     );

    std::cout<<"max_temperature="<<max_temperature<<"\n";
    std::cout<<"energy_density ="<<energy_density <<"\n";
    std::cout<<"power_density  ="<<power_density  <<"\n";
    std::cout<<"efficiency     ="<<efficiency     <<"\n";

    return 0;
}

#include "common.cc"
