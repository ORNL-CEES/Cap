#define BOOST_TEST_MODULE TestLissajousCurve
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/utils.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <iterator>
#include <algorithm>

namespace cap {

void print_headers(std::string const & mode, std::ostream & os)
{
    os<<"# lissajous curve (mode="<<mode<<")\n";
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
get_evolve_one_time_step(std::string const & mode, std::shared_ptr<boost::property_tree::ptree const> database)
{
    std::vector<int   > const harmonics  = cap::to_vector<int   >(database->get<std::string>("harmonics" ));
    std::vector<double> const amplitudes = cap::to_vector<double>(database->get<std::string>("amplitudes"));
    std::vector<double> const phases     = cap::to_vector<double>(database->get<std::string>("phases"    ));
    double const frequency   = database->get<double>("frequency");
    double const pi          = std::acos(-1.0);
    int    const n_harmonics = harmonics.size();
    std::vector<double> angular_frequencies(n_harmonics);
    std::transform(
        harmonics.begin(), harmonics.end(),
        angular_frequencies.begin(),
        [&pi,&frequency](int k) { return 2.0 * pi * k * frequency; }
        );
    auto compute_ac_excitation_signal =
        [n_harmonics,amplitudes,angular_frequencies,phases](double time)
        {
            double excitation_signal = 0.0;
            for (int n = 0; n < n_harmonics; ++n)
                excitation_signal += amplitudes[n] * std::sin(angular_frequencies[n] * time + phases[n]);
            return excitation_signal;
        };

    if (mode.compare("galvanostatic") == 0) {
        double const dc_current = database->get<double>("dc_current");
        return [dc_current,compute_ac_excitation_signal](double const time, double const time_step, std::shared_ptr<cap::EnergyStorageDevice> dev)
            { dev->evolve_one_time_step_changing_current(time_step, dc_current+compute_ac_excitation_signal(time)); };
    } else if (mode.compare("potentiostatic") == 0) {
        double const dc_voltage = database->get<double>("dc_voltage");
        return [dc_voltage,compute_ac_excitation_signal](double const time, double const time_step, std::shared_ptr<cap::EnergyStorageDevice> dev)
            { dev->evolve_one_time_step_changing_voltage(time_step, dc_voltage+compute_ac_excitation_signal(time)); };
    } else {
        throw std::runtime_error("invalid EIS measurement mode "+mode);
    }
}

void measure_lissajous_curve(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    std::string const mode = database->get<std::string>("mode");
    auto evolve_one_time_step = get_evolve_one_time_step(mode, database);

    int    const cycles          = database->get<int   >("cycles"         );
    int    const ignore_cycles   = database->get<int   >("ignore_cycles"  );
    int    const steps_per_cycle = database->get<int   >("steps_per_cycle");
    double const frequency       = database->get<double>("frequency"      );
    double const initial_voltage = database->get<double>("initial_voltage");

    double       time      = 0.0;
    double const time_step = 1.0 / frequency / steps_per_cycle;
    dev->reset_voltage(initial_voltage);
    for (int n = 0; n < cycles*steps_per_cycle; ++n)
    {
          time += time_step;
          evolve_one_time_step(time, time_step, dev);
          if (n >= ignore_cycles*steps_per_cycle)
              report(time, dev, os);
    }

}



std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
recycle(std::istream & is)
{
//    char buffer[256];
//    is.getline(buffer, 256);
//    std::cout<<buffer;
    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> voltage;
    while (is.good())
    {
        time.emplace_back(std::nan(""));
        is>>time.back();
        if (!is.eof()) {
            current.emplace_back(std::nan(""));
            voltage.emplace_back(std::nan(""));
            is>>current.back()>>voltage.back();
            std::cout<<boost::format("  time=%22.15e  current=%22.15e  voltage=%22.15e  \n")
                % time.back()
                % current.back()
                % voltage.back()
                ;
        } else {
            time.pop_back();
        }
    }
    std::cout<<"size="<<time.size()<<"\n";
    std::cout<<"total time="<<time.back()-time.front()<<"\n";
    return std::make_tuple(time, current, voltage);
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_lissajous_curve )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_lissajous_curve", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // measure its impedance
    std::fstream fout;
    fout.open("lissajous_curve_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> lissajous_curve_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("lissajous_curve"));
    cap::measure_lissajous_curve(device, lissajous_curve_database, fout);

    std::stringstream ss;
    cap::measure_lissajous_curve(device, lissajous_curve_database, ss);
    std::vector<double> current;
    std::vector<double> voltage;
    std::tie(std::ignore, current, voltage) = cap::recycle(ss);

    fout.close();
}    
