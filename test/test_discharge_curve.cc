#define BOOST_TEST_MODULE TestDischargeCurve
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/utils.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_fft_real.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <iterator>
#include <algorithm>



namespace cap {

void print_headers(std::string const & discharge_mode, std::ostream & os)
{
    os<<"# discharge curve (discharge_mode="<<discharge_mode<<")\n";
    os<<boost::format( "# %22s  %22s  %22s  \n")
        % "time t [seconds]"
        % "current I [amperes]"
        % "voltage U [volts]"
        ;
}



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
}



std::function<void(double const, double const, std::shared_ptr<cap::EnergyStorageDevice>)>
get_evolve_one_time_step(std::string const & discharge_mode, std::shared_ptr<boost::property_tree::ptree const> database)
{
    if (discharge_mode.compare("constant_current") == 0) {
        double const discharge_current = database->get<double>("discharge_current");
        return [discharge_current](double const, double const time_step, std::shared_ptr<cap::EnergyStorageDevice> dev)
            { dev->evolve_one_time_step_constant_current(time_step, -discharge_current); };
    } else if (discharge_mode.compare("constant_power") == 0) {
        double const discharge_power = database->get<double>("discharge_power");
        return [discharge_power](double const, double const time_step, std::shared_ptr<cap::EnergyStorageDevice> dev)
            { dev->evolve_one_time_step_constant_power(time_step, -discharge_power); };
    } else if (discharge_mode.compare("constant_load") == 0) {
        double const discharge_load = database->get<double>("discharge_load");
        return [discharge_load](double const, double const time_step, std::shared_ptr<cap::EnergyStorageDevice> dev)
            { dev->evolve_one_time_step_constant_load(time_step, discharge_load); };
    } else {
        throw std::runtime_error("invalid discharge mode "+discharge_mode);
    }
}



std::function<bool(double const, std::shared_ptr<cap::EnergyStorageDevice>)>
get_end_criterion(std::string const & discharge_stop_at, std::shared_ptr<boost::property_tree::ptree const> database)
{
    if (discharge_stop_at.compare("voltage_lower_than") == 0) {
        double const discharge_voltage_limit = database->get<double>("discharge_voltage_limit");
        return [discharge_voltage_limit](double const, std::shared_ptr<cap::EnergyStorageDevice> dev)
            { double voltage; dev->get_voltage(voltage); return (voltage <= discharge_voltage_limit); };
    } else if (discharge_stop_at.compare("time_greater_than") == 0) {
        double const discharge_duration = database->get<double>("discharge_duration");
        return [discharge_duration](double const time, std::shared_ptr<cap::EnergyStorageDevice>)
            { return (time >= discharge_duration); };
    } else {
        throw std::runtime_error("invalid test stop at "+discharge_stop_at);
    }
}



void discharge_curve(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    std::string const discharge_mode    = database->get<std::string>("discharge_mode"   );
    auto evolve_one_time_step = get_evolve_one_time_step(discharge_mode, database);
    std::string const discharge_stop_at = database->get<std::string>("discharge_stop_at");
    auto end_criterion = get_end_criterion(discharge_stop_at, database);

    double const max_discharge_time = database->get<double>("max_discharge_time");
    double const time_step          = database->get<double>("time_step"         );
    double const initial_voltage    = database->get<double>("initial_voltage"   );
    dev->reset_voltage(initial_voltage);

    double time = 0.0;
    print_headers(discharge_mode, os);
    report(time, dev, os);
    while (time < max_discharge_time)
    {
        evolve_one_time_step(time, time_step, dev);
        time += time_step;
        report(time, dev, os);
        if (end_criterion(time, dev))
            break;
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_discharge_curve )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_discharge_curve", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // measure discharge curve
    std::fstream fout;
    fout.open("discharge_curve_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> discharge_curve_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("discharge_curve"));
    cap::discharge_curve(device, discharge_curve_database, fout);

}    
