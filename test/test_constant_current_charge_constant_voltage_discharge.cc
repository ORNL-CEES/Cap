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
    std::ignore = argc;
    std::ignore = argv;
    cap::SuperCapacitorProblem<2> super_capacitor(initialize_database());

    // SETTING PROBLEM PARAMETERS
    std::shared_ptr<boost::property_tree::ptree> in(new boost::property_tree::ptree);
    put_default_parameters(*in);
    in->put("test_case", 1);
    in->put("time_step",       "0.5");
    in->put("initial_time",    "0.0");
    in->put("final_time",    "900.0");
    in->put("boundary_values.charge_current_density", 324.65);
    in->put("boundary_values.discharge_potential",      1.4 );

    // SOLVING THE PROBLEM
    std::shared_ptr<boost::property_tree::ptree> out(new boost::property_tree::ptree);
    super_capacitor.run(in, out);

    // POSTPROCESSING QUANTITIES OF INTEREST
    std::vector<double> max_temperature = 
        cap::to_vector<double>(out->get<std::string>("max_temperature"));

    return 0;
}

#include "common.cc"
