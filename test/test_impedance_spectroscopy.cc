#define BOOST_TEST_MODULE TestImpedanceSpectroscopy
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

std::vector<double>
compute_signal_to_noise_ratio2(std::vector<double> const & data, std::vector<int> const & components)
{
    size_t const n = data.size();
    std::vector<size_t> excited_harmonics;
    std::vector<size_t> unexcited_harmonics;
    for (size_t i = 1; i < n/2; ++i)
        if (std::find(components.begin(), components.end(), i) != components.end())
            excited_harmonics.push_back(i);
        else
            unexcited_harmonics.push_back(i);
    if (components.size() != excited_harmonics.size())
        throw std::runtime_error("check your harmonics");
    if (excited_harmonics.size() + unexcited_harmonics.size() != n/2 - 1)
        throw std::runtime_error("il en manque une");
    
    double dc_power = std::pow(data[0], 2);
    double ac_power = 0.0;
    for (size_t k : excited_harmonics)
        ac_power += std::pow(data[k], 2) + std::pow(data[n-k], 2);
    double noise_power = std::pow(data[n/2], 2);
    for (size_t k : unexcited_harmonics)
        noise_power += std::pow(data[k], 2) + std::pow(data[n-k], 2);
    double total_power = dc_power + ac_power + noise_power;
    std::ignore = total_power;
    // TODO:  definition is a bit shaky for multiple frequency
    // should compute for each individual excited freq vs nearby unex freq
    std::vector<double> ratios;
    double harmonics_power;
    for (size_t k : excited_harmonics) {
        harmonics_power = std::pow(data[k], 2) + std::pow(data[n-k], 2);
        ratios.push_back(harmonics_power / noise_power);
    }
    return ratios;
}



double compute_signal_to_noise_ratio(std::vector<double> const & data, size_t i)
{
    return compute_signal_to_noise_ratio2(data, std::vector<int>{static_cast<int>(i)})[0];
}



std::vector<std::complex<double>>
measure_impedance2(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database)
{
    std::string const mode = database->get<std::string>("mode");
    if ((mode.compare("galvanostatic") != 0)
        && (mode.compare("potentiostatic") != 0))
        throw std::runtime_error("invalide mode "+mode+" in EIS measurement");

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
        [n_harmonics,&amplitudes,&angular_frequencies,&phases](double time)
        {
            double excitation_signal = 0.0;
            for (int n = 0; n < n_harmonics; ++n)
                excitation_signal += amplitudes[n] * std::sin(angular_frequencies[n] * time + phases[n]);
            return excitation_signal;
        };
    auto evolve_one_time_step =
        (
            (mode.compare("galvanostatic") == 0)
        ?
            [](std::shared_ptr<cap::EnergyStorageDevice> device, double time, double time_step, double & voltage, double & current, decltype(compute_ac_excitation_signal) & compute_ac_excitation_signal)
            {
                current += compute_ac_excitation_signal(time);
                device->evolve_one_time_step_changing_current(time_step, current);
                device->get_voltage(voltage);
            }
        :
            [](std::shared_ptr<cap::EnergyStorageDevice> device, double time, double time_step, double & voltage, double & current, decltype(compute_ac_excitation_signal) & compute_ac_excitation_signal)
            {
                voltage += compute_ac_excitation_signal(time);
                device->evolve_one_time_step_changing_voltage(time_step, voltage);
                device->get_current(current);
            }
        );

    double const dc_voltage      = database->get<double>("dc_voltage"     );
    int    const cycles          = database->get<int   >("cycles"         );
    int    const ignore_cycles   = database->get<int   >("ignore_cycles"  );
    int    const steps_per_cycle = database->get<int   >("steps_per_cycle");
    double const initial_voltage = dc_voltage + compute_ac_excitation_signal(0.0);
    std::vector<int> const powers_of_two =
        { 1, 2 ,4, 8, 16,
          32, 64, 128, 256, 512,
          1024, 2048, 4096, 8192, 16384,
          32768, 65536, 131072, 262144, 524288
        }; // first 20 powers of two
    int const exponent =
        std::distance(
            powers_of_two.begin(),
            std::find(powers_of_two.begin(), powers_of_two.end(),
                (cycles - ignore_cycles) * steps_per_cycle)
            );
    if (exponent == static_cast<int>(powers_of_two.size()))
        throw std::runtime_error("(cycles-ignore_cycles)*steps_per_cycle must be a power of two");
    if (powers_of_two[exponent] != (1 << exponent))
        throw std::runtime_error("something wrong with the powers of two");

    double       time      = 0.0;
    double const time_step = 1.0 / frequency / steps_per_cycle;
    dev->reset_voltage(initial_voltage);
    std::vector<double> voltage;
    std::vector<double> current;
    voltage.reserve(cycles*steps_per_cycle);
    current.reserve(cycles*steps_per_cycle);
    for (int n = 0; n < cycles*steps_per_cycle; ++n)
    {
          if (n == ignore_cycles*steps_per_cycle) {
              voltage.clear();
              current.clear();
          }
          voltage.push_back(dc_voltage);
          current.push_back(0.0);
          time += time_step;
          evolve_one_time_step(dev, time, time_step, voltage.back(), current.back(), compute_ac_excitation_signal);
    }
    size_t const n = voltage.size();
    gsl_fft_real_radix2_transform(&(voltage[0]), 1, n);
    gsl_fft_real_radix2_transform(&(current[0]), 1, n);

    // compute signal-to-noise ratio
    std::vector<int> dummy_components(harmonics);
    for(int & x : dummy_components)
         x *= (cycles-ignore_cycles);

    std::cout<<"voltage signal-to-noise ratio = "<<compute_signal_to_noise_ratio2(voltage, dummy_components)[0]<<"\n";
    std::cout<<"current signal-to-noise ratio = "<<compute_signal_to_noise_ratio2(current, dummy_components)[0]<<"\n";

    std::vector<std::complex<double>> results;
    for (auto const & k : dummy_components)
        results.push_back(
            std::complex<double>(voltage[k], voltage[n-k])
            /
            std::complex<double>(current[k], current[n-k])
        );
    return results;
}



std::complex<double>
measure_impedance(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database)
{
    double const frequency       = database->get<double>("frequency"      );
    int    const cycles          = database->get<int   >("cycles"         );
    int    const ignore_cycles   = database->get<int   >("ignore_cycles"  );
    int    const steps_per_cycle = database->get<int   >("steps_per_cycle");
    double const dc_voltage      = database->get<double>("dc_voltage"     );
    double const ac_component    = database->get<double>("ac_component"   );
    double const initial_voltage = dc_voltage;
    std::vector<int> const powers_of_two =
        { 1, 2 ,4, 8, 16, 
          32, 64, 128, 256, 512, 
          1024, 2048, 4096, 8192, 16384, 
          32768, 65536, 131072, 262144, 524288
        }; // first 20 powers of two
    if (find(powers_of_two.begin(), powers_of_two.end(), (cycles-ignore_cycles)*steps_per_cycle) == powers_of_two.end())
        throw std::runtime_error("(cycles-1)*steps_per_cycle must be a power of two");

    double       time      = 0.0;
    double const time_step = 1.0 / frequency / steps_per_cycle;
    double const pi        = std::acos(-1.0);
    double const angular_frequency = 2.0 * pi * frequency;
    dev->reset_voltage(initial_voltage);
    double voltage;
    double current;
    std::vector<double> excitation;
    std::vector<double> response;
    for (int n = 0; n < cycles*steps_per_cycle; ++n)
    {
          time += time_step;
          voltage = dc_voltage + ac_component * std::sin(angular_frequency * time);
          dev->evolve_one_time_step_changing_voltage(time_step, voltage);
          dev->get_current(current);
          if (n >= ignore_cycles*steps_per_cycle)
          {
              excitation.push_back(voltage);
              response  .push_back(current);
          }
    }
    gsl_fft_real_radix2_transform(&(excitation[0]), 1, excitation.size());
    gsl_fft_real_radix2_transform(&(response  [0]), 1, response  .size());
//    std::cout<<"excitation signal-to-noise ratio = "<<compute_signal_to_noise_ratio(excitation, cycles-ignore_cycles)<<"\n";
//    std::cout<<"response signal-to-noise ratio = "<<compute_signal_to_noise_ratio(response, cycles-ignore_cycles)<<"\n";

    std::complex<double> const impedance =
        std::complex<double>(excitation[(cycles-ignore_cycles)], excitation[excitation.size()-(cycles-ignore_cycles)])
        /
        std::complex<double>(response  [(cycles-ignore_cycles)], response  [response  .size()-(cycles-ignore_cycles)]);
    return impedance;
}



void impedance_spectroscopy(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const frequency_lower_limit = database->get<double>("frequency_lower_limit");
    double const frequency_upper_limit = database->get<double>("frequency_upper_limit");
    int    const steps_per_decade      = database->get<int   >("steps_per_decade"     );
    double const pi                    = boost::math::constants::pi<double>();
    std::shared_ptr<boost::property_tree::ptree> tmp =
        std::make_shared<boost::property_tree::ptree>(*database);
    for (double frequency = frequency_lower_limit; frequency <= frequency_upper_limit; frequency *= std::pow(10.0, 1.0/steps_per_decade))
    {
        tmp->put("frequency", frequency);
        std::complex<double> impedance =
            measure_impedance(dev, tmp);
        os<<boost::format( "  %20.15e  %20.15e  %20.15e  %20.15e  %20.15e  \n")
            % frequency
            % impedance.real()
            % impedance.imag()
            % std::abs(impedance)
            % (std::arg(impedance) * 180.0 / pi)
            ;
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_measure_impedance2 )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_impedance_spectroscopy", *input_database);

    std::string const type = input_database->get<std::string>("device.type");
    // current test will only work for series or parallel rc circuit
    if ((type.compare("SeriesRC") != 0) && (type.compare("ParallelRC") != 0))
        throw std::runtime_error("test measure impedance check not implemented for "+type);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    double const frequency_lower_limit = input_database->get<double>("impedance_spectroscopy.frequency_lower_limit");
    double const frequency_upper_limit = input_database->get<double>("impedance_spectroscopy.frequency_upper_limit");
    std::shared_ptr<boost::property_tree::ptree> impedance_spectroscopy_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("impedance_spectroscopy"));

    std::cout<<"###############\n";
    std::cout<<"# MULTI FREQ  #\n";
    std::cout<<"###############\n";
    double const base_frequency = 0.5 * (frequency_lower_limit + frequency_upper_limit);
    impedance_spectroscopy_database->put("frequency", base_frequency);
    for (auto const & impedance : cap::measure_impedance2(device, impedance_spectroscopy_database))
        std::cout<<impedance<<"\n";

    std::cout<<"###############\n";
    std::cout<<"# SINGLE FREQ #\n";
    std::cout<<"###############\n";
    std::vector<int> const harmonics = cap::to_vector<int>(input_database->get<std::string>("impedance_spectroscopy.harmonics"));
    for (int k : harmonics) {
        impedance_spectroscopy_database->put("frequency", k*base_frequency);
        std::cout<<cap::measure_impedance(device, impedance_spectroscopy_database)<<"\n";
    }
}




BOOST_AUTO_TEST_CASE( test_measure_impedance )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_impedance_spectroscopy", *input_database);

    std::string const type = input_database->get<std::string>("device.type");
    // current test will only work for series or parallel rc circuit
    if ((type.compare("SeriesRC") != 0) && (type.compare("ParallelRC") != 0))
        throw std::runtime_error("test measure impedance check not implemented for "+type);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    input_database->put("impedance_spectroscopy.harmonics", "1");

    double const frequency_lower_limit = input_database->get<double>("impedance_spectroscopy.frequency_lower_limit");
    double const frequency_upper_limit = input_database->get<double>("impedance_spectroscopy.frequency_upper_limit");
    int    const steps_per_decade      = input_database->get<int   >("impedance_spectroscopy.steps_per_decade"     );
    double const percent_tolerance     = input_database->get<double>("impedance_spectroscopy.percent_tolerance"    );
    double const series_resistance     = input_database->get<double>("device.series_resistance"                    );
    double const parallel_resistance   = input_database->get<double>("device.parallel_resistance"                  );
    double const capacitance           = input_database->get<double>("device.capacitance"                          );
    double const pi                    = boost::math::constants::pi<double>();

    std::shared_ptr<boost::property_tree::ptree> impedance_spectroscopy_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("impedance_spectroscopy"));
    std::fstream fout("computed_vs_exact_impedance_spectroscopy_data", std::fstream::out);
    for (double frequency = frequency_lower_limit; frequency < frequency_upper_limit; frequency *= std::pow(10.0, 1.0/steps_per_decade))
    {
        double const angular_frequency = 2.0 * pi * frequency;
        impedance_spectroscopy_database->put("frequency", frequency);
        std::complex<double> computed_impedance = cap::measure_impedance2(device, impedance_spectroscopy_database)[0];
        std::complex<double> exact_impedance    =
            (
                (type.compare("SeriesRC") == 0)
                ?
                series_resistance + 1.0 / std::complex<double>(0.0, capacitance * angular_frequency)
                :
                series_resistance + parallel_resistance / std::complex<double>(1.0, parallel_resistance * capacitance * angular_frequency)
            );
        fout<<boost::format("  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  \n")
            % frequency
            % computed_impedance.real()
            % computed_impedance.imag()
            % std::abs(computed_impedance)
            % (std::arg(computed_impedance) * 180.0 / pi)
            % exact_impedance.real()
            % exact_impedance.imag()
            % std::abs(exact_impedance)
            % (std::arg(exact_impedance) * 180.0 / pi)
            ;
//          BOOST_CHECK_CLOSE(computed_impedance.real(), exact_impedance.real(), percent_tolerance);
//          BOOST_CHECK_CLOSE(computed_impedance.imag(), exact_impedance.imag(), percent_tolerance);
          BOOST_CHECK_CLOSE(std::abs(computed_impedance), std::abs(exact_impedance), percent_tolerance);
          BOOST_CHECK_CLOSE(std::arg(computed_impedance), std::arg(exact_impedance), percent_tolerance);
    }

}



BOOST_AUTO_TEST_CASE( test_impedance_spectroscopy )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_impedance_spectroscopy", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // measure its impedance
    std::fstream fout;
    fout.open("impedance_spectroscopy_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> impedance_spectroscopy_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("impedance_spectroscopy"));
    cap::impedance_spectroscopy(device, impedance_spectroscopy_database, fout);

}    
