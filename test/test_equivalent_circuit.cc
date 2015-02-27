#define BOOST_TEST_MODULE TestEquivalentCircuit
#define BOOST_TEST_MAIN
#include <cap/super_capacitor.h>
#include <cap/equivalent_circuit.h>
#include <cap/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <functional>
#include <algorithm>
#include <memory>


void report_time_current_voltage(std::vector<double> const & time,
    std::vector<double> const & current,
    std::vector<double> const & voltage,
    std::ostream & os = std::cout)
{
    std::size_t const n = time.size();
    if ((current.size() != n)
        || (voltage.size() != n))
        throw std::runtime_error("vector size mismatch");
    for (std::size_t i = 0; i < n; ++i)
        os<<boost::format("  %10.1f  %10.7f  %10.5f  \n")
            % time[i]
            % current[i]
            % voltage[i]
            ;
}

BOOST_AUTO_TEST_CASE( test_report_time_current_voltage )
{
    std::vector<double> t(1);
    std::vector<double> i(2);
    std::vector<double> u(1);
    BOOST_CHECK_THROW( report_time_current_voltage(t, i, u), std::runtime_error );
}

// forward declaration (implementation is in "common.cc")
std::shared_ptr<boost::property_tree::ptree> initialize_database();
void put_default_parameters(boost::property_tree::ptree & params);

BOOST_AUTO_TEST_CASE( test_equivalent_circuit )
{
    std::shared_ptr<boost::property_tree::ptree> in(new boost::property_tree::ptree);
    put_default_parameters(*in);
    in->put("test_case",       2  );
    in->put("time_step",       0.1);
    in->put("initial_time",    0.0);
    in->put("final_time",     15.0);
    in->put("max_cycles",      1  );
    in->put("debug.solution_fields", "solid_potential,liquid_potential,overpotential,joule_heating");
    in->put("debug.solution_fluxes", "liquid_current_density,solid_current_density");
    in->put("debug.material_properties", "specific_capacitance,liquid_electrical_conductivity");
    in->put("debug.material_ids", true);

    in->put("material_properties.electrode_void_volume_fraction", 0.55);

    std::shared_ptr<boost::property_tree::ptree> out(new boost::property_tree::ptree);
    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> voltage;
    std::vector<double> power;
    std::vector<double> energy;
    std::vector<double> heat_production;
    std::vector<double> thermal_losses;
    std::vector<std::string> capacitor_state;

    std::cout<<"########  BEGIN HIGH FIDELITY MODEL  ########\n";
    std::shared_ptr<boost::property_tree::ptree> database = initialize_database();
    database->put("verbose", false);
    cap::SuperCapacitorProblem<2> high_fidelity_model(database);
    high_fidelity_model.run(in, out);
    time    = cap::to_vector<double>(out->get<std::string>("time")   );
    current = cap::to_vector<double>(out->get<std::string>("current"));
    voltage = cap::to_vector<double>(out->get<std::string>("voltage"));
    capacitor_state = cap::to_vector<std::string>(out->get<std::string>("capacitor_state"));
    heat_production = cap::to_vector<double>(out->get<std::string>("heat_production"));
    power = cap::to_vector<double>(out->get<std::string>("power_density"));
    energy = cap::to_vector<double>(out->get<std::string>("energy_density"));
    thermal_losses = cap::to_vector<double>(out->get<std::string>("thermal_energy_loss"));
    report_time_current_voltage(time, current, voltage);
    std::cout<<"########  END HIGH FIDELITY MODEL  ########\n";

    // clear it entirely to avoid silly mistakes
    double const  mass = out->get<double>("mass"); // need to get the mass from here for now
    out->clear();

    std::cout<<"########  BEGIN LOW FIDELITY MODEL  ########\n";
    cap::EquivalentCircuit low_fidelity_model(database);
    low_fidelity_model.run(in, out);
    time    = cap::to_vector<double>(out->get<std::string>("time")   );
    current = cap::to_vector<double>(out->get<std::string>("current"));
    voltage = cap::to_vector<double>(out->get<std::string>("voltage"));
    capacitor_state = cap::to_vector<std::string>(out->get<std::string>("capacitor_state"));
    heat_production = cap::to_vector<double>(out->get<std::string>("heat_production"));
    // compute these by hands for now
    std::size_t const n = time.size();
    power.resize(n);
    energy.resize(n);
    thermal_losses.resize(n);
    std::transform(voltage.begin(), voltage.end(), current.begin(), power.begin(), 
        [&mass](double const & U, double const & I) { return U * I / mass; });
    cap::compute_energy(capacitor_state, time, power, energy);
    cap::compute_thermal_energy_losses(capacitor_state, time,
        heat_production, thermal_losses);

    report_time_current_voltage(time, current, voltage);
    std::cout<<"########  END LOW FIDELITY MODEL  ########\n";


}

#include "common.cc"
