#include <cap/super_capacitor.h>
#include <cap/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>

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
    in->put("time_step",       5.0);
    in->put("initial_time",    0.0);
    in->put("final_time",    600.0);
    in->put("max_cycles",    100  );
    in->put("boundary_values.charge_current_density",     324.65);      
    in->put("boundary_values.discharge_current_density", -324.65);      

    // SOLVING THE PROBLEM
    std::shared_ptr<boost::property_tree::ptree> out(new boost::property_tree::ptree);
    super_capacitor.run(in, out);

    // POSTPROCESSING QUANTITIES OF INTEREST
    std::vector<double> max_temperature = cap::to_vector<double>(out->get<std::string>("max_temperature"));
    std::vector<double> heat_production = cap::to_vector<double>(out->get<std::string>("heat_production"));
    std::vector<double> current         = cap::to_vector<double>(out->get<std::string>("current")        );
    std::vector<double> voltage         = cap::to_vector<double>(out->get<std::string>("voltage")        );
    std::vector<double> time            = cap::to_vector<double>(out->get<std::string>("time")           );
    double volume = out->get<double>("volume");
    double mass   = out->get<double>("mass");
    std::vector<std::string> capacitor_state = cap::to_vector<std::string>(out->get<std::string>("capacitor_state"));
    std::vector<int>         cycle           = cap::to_vector<int>        (out->get<std::string>("cycle")          );

    std::size_t n = time.size();
    std::vector<double> power(n);
    std::transform(voltage.begin(), voltage.end(), current.begin(), power.begin(), std::multiplies<double>());
    std::vector<double> energy(n);
    energy[0] = 0.0;
    for (std::size_t i = 1; i < n; ++i)
        if (capacitor_state[i].compare(capacitor_state[i-1]) == 0)
            energy[i] = energy[i-1] + 0.5 * (time[i] - time[i-1]) * (power[i] + power[i-1]);
        else
            energy[i] = energy[i-1];
    std::vector<double> efficiency(n);
    std::transform(power.begin(), power.end(), heat_production.begin(), efficiency.begin(), 
        [](double const & P, double const & Q) { return 100.0 * (std::abs(P) - Q) / std::abs(P); });



    for (std::size_t i = 0; i < n; ++i) {
        std::cout<<i<<"  "
            <<time[i]<<"  "
            <<capacitor_state[i]<<"  "
            <<cycle[i]<<"  "
            <<current[i]<<"  "
            <<voltage[i]<<"  "
            <<max_temperature[i]<<"  "
            <<heat_production[i]<<"  "
            <<power[i]<<"  "
            <<efficiency[i]<<"  "
            <<energy[i]<<"  "
            <<"\n";
    } // end for i
    std::cout
        <<volume<<"  "
        <<mass<<"\n";


    auto minmax_energy = std::minmax_element(energy.begin(), energy.end());
    double seconds_per_hour = 3600.0;
    double energy_density = (*minmax_energy.second - *minmax_energy.first) / (mass * seconds_per_hour);
    std::cout<<energy_density<<"\n";

    double power_density = std::accumulate(power.begin(), power.end(), 0.0, 
        [](double const & sum, double const & val) { return sum + std::abs(val); }) / (static_cast<double>(power.size()) * mass);
    std::cout<<power_density<<"\n";

    return 0;
}

#include "common.cc"
