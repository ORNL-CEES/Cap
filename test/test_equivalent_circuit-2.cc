#define BOOST_TEST_MODULE TestEquivalentCircuit
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/equivalent_circuit.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>

namespace cap {

void report(double const time, std::shared_ptr<EnergyStorageDevice const> dev, std::ostream & os)
{
    double current;
    dev->get_current(current);
    double voltage;
    dev->get_voltage(voltage);
    os<<boost::format("  %22.15e  %22.15e  %22.15e  \n")
        % time
        % current
        % voltage
        ;
    // debug output of the solution fields and fluxes
    std::ostream null_sink(0);
    dev->print_data(null_sink);
}



void cycling_charge_discharge(std::shared_ptr<boost::property_tree::ptree const> database, std::shared_ptr<cap::EnergyStorageDevice> dev, std::ostream & os)
{
    double const time_step           = database->get<double>("time_step"               );
    double const initial_voltage     = database->get<double>("initial_voltage"         );
    double const charge_current      = database->get<double>("charge_current"          );
    double const discharge_current   = database->get<double>("discharge_current"       );
    double const voltage_upper_limit = database->get<double>("voltage_upper_limit"     );
    double const voltage_lower_limit = database->get<double>("voltage_lower_limit"     );
    double const charge_hold_time    = database->get<double>("charge_hold_time"   , 0.0);
    double const discharge_hold_time = database->get<double>("discharge_hold_time", 0.0);
    double const charge_rest_time    = database->get<double>("charge_rest_time"   , 0.0);
    double const discharge_rest_time = database->get<double>("discharge_rest_time", 0.0);
    int    const cycles              = database->get<int   >("cycles"                  );
    dev->reset_voltage(initial_voltage);
    double voltage = initial_voltage;
    double time = 0.0;
    double tick;
    for (int n = 0; n < cycles; ++n)
    {
        // glavanostatic charge
        while (voltage < voltage_upper_limit)
        {
            dev->evolve_one_time_step_constant_current(time_step, charge_current);
            dev->get_voltage(voltage);
            time += time_step;
            report(time, dev, os);
        }
        // potentiostatic hold
        tick = time;
        while (time - tick < charge_hold_time)
        {
            dev->evolve_one_time_step_constant_voltage(time_step, voltage);
            time += time_step;
            report(time, dev, os);
        }
        // rest at open circuit
        tick = time;
        while (time - tick < charge_rest_time)
        {
            dev->evolve_one_time_step_constant_current(time_step, 0.0);
            time += time_step;
            report(time, dev, os);
        }
        // glavanostatic discharge
        while (voltage > voltage_lower_limit)
        {
            dev->evolve_one_time_step_constant_current(time_step, discharge_current);
            dev->get_voltage(voltage);
            time += time_step;
            report(time, dev, os);
        }
        // potentiostatic hold
        tick = time;
        while (time - tick < discharge_hold_time)
        {
            dev->evolve_one_time_step_constant_voltage(time_step, voltage);
            time += time_step;
            report(time, dev, os);
        }
        // rest at open circuit
        tick = time;
        while (time - tick < discharge_rest_time)
        {
            dev->evolve_one_time_step_constant_current(time_step, 0.0);
            time += time_step;
            report(time, dev, os);
        }

    }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_equivalent_circuit )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_equivalent_circuit", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // 
    std::shared_ptr<boost::property_tree::ptree> not_empty_database =
        std::make_shared<boost::property_tree::ptree>();
    not_empty_database->put("something", "not_empty");
    BOOST_CHECK_THROW(cap::compute_equivalent_circuit(device_database, not_empty_database), std::runtime_error);

    // build an equivalent resistor capacitor circuit
    std::shared_ptr<boost::property_tree::ptree> equivalent_circuit_database =
        std::make_shared<boost::property_tree::ptree>();
    cap::compute_equivalent_circuit(device_database, equivalent_circuit_database);
    std::shared_ptr<cap::EnergyStorageDevice> equivalent_circuit =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(equivalent_circuit_database));
    
    std::shared_ptr<boost::property_tree::ptree> cycling_charge_discharge_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("cycling_charge_discharge"));

    std::fstream fout;
    fout.open("equivalent_circuit_data_lf", std::fstream::out);
    cap::cycling_charge_discharge(cycling_charge_discharge_database, equivalent_circuit, fout);

    fout.close();
    fout.open("equivalent_circuit_data_hf", std::fstream::out);
    cap::cycling_charge_discharge(cycling_charge_discharge_database, device            , fout);
    
}
