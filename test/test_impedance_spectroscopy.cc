#define BOOST_TEST_MODULE TestImpedanceSpectroscopy
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
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



namespace cap {

std::complex<double>
measure_impedance(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout);

void bar(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const frequency_lower_limit = database->get<double>("frequency_lower_limit");
    double const frequency_upper_limit = database->get<double>("frequency_upper_limit");
    double const ratio                 = database->get<double>("ratio"                );
    double const pi                    = boost::math::constants::pi<double>();
    std::shared_ptr<boost::property_tree::ptree> tmp =
        std::make_shared<boost::property_tree::ptree>(*database);
    for (double frequency = frequency_lower_limit; frequency <= frequency_upper_limit; frequency *= ratio)
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



std::complex<double>
measure_impedance(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os)
{
    double const frequency       = database->get<double>("frequency"      );
    double const amplitude       = database->get<double>("amplitude"      );
    int    const cycles          = database->get<int   >("cycles"         );
    int    const steps_per_cycle = database->get<int   >("steps_per_cycle");
    double const initial_voltage = 0.0;
    std::vector<int> const powers_of_two =
        { 1, 2 ,4, 8, 16, 
          32, 64, 128, 256, 512, 
          1024, 2048, 4096, 8192, 16384, 
          32768, 65536, 131072, 262144, 524288
        }; // first 20 powers of two
    if (find(powers_of_two.begin(), powers_of_two.end(), (cycles-1)*steps_per_cycle) == powers_of_two.end())
        throw std::runtime_error("(cycles-1)*steps_per_cycle must be a power of two");

    double       time      = 0.0;
    double const time_step = 1.0 / frequency / steps_per_cycle;
    double const phase     = std::asin(initial_voltage / amplitude);
    double const pi        = std::acos(-1.0);
    dev->reset_voltage(initial_voltage);
    double voltage;
    double current;
    std::vector<double> excitation;
    std::vector<double> response;
    for (int n = 0; n < cycles*steps_per_cycle; ++n)
    {
          time += time_step;
          voltage = amplitude*std::sin(2.0*pi*frequency*time+phase);
          dev->evolve_one_time_step_constant_voltage(time_step, voltage);
          dev->get_current(current);
          if (n >= steps_per_cycle)
          {
              excitation.push_back(voltage);
              response  .push_back(current);
          }
//          os<<boost::format("  %20.15e  %20.15e  %20.15e  \n")
//              % time
//              % current
//              % voltage
//              ;
    }
    gsl_fft_real_radix2_transform(&(excitation[0]), 1, excitation.size());
    gsl_fft_real_radix2_transform(&(response  [0]), 1, response  .size());
    std::complex<double> const impedance =
        std::complex<double>(excitation[(cycles-1)], excitation[excitation.size()-(cycles-1)])
        /
        std::complex<double>(response  [(cycles-1)], response  [response  .size()-(cycles-1)]);
    return impedance;
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_measure_impedance )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_impedance_spectroscopy", *input_database);

    std::string const type = input_database->get<std::string>("device.type");
    if ((type.compare("SeriesRC") != 0) && (type.compare("ParallelRC") != 0))
        input_database->put("device.type", "ParallelRC");
    std::string const duh = input_database->get<std::string>("device.type");
    if (duh.compare("NoName") == 0)
        throw std::runtime_error("duh");
    

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    double const frequency_lower_limit = input_database->get<double>("impedance_spectroscopy.frequency_lower_limit");
    double const frequency_upper_limit = input_database->get<double>("impedance_spectroscopy.frequency_upper_limit");
    double const ratio                 = input_database->get<double>("impedance_spectroscopy.ratio"                );
    double const series_resistance     = input_database->get<double>("device.series_resistance"  );
    double const parallel_resistance   = input_database->get<double>("device.parallel_resistance");
    double const capacitance           = input_database->get<double>("device.capacitance"        );
    double const pi                    = boost::math::constants::pi<double>();
    double const percent_tolerance     = 1.0e-1;
    std::shared_ptr<boost::property_tree::ptree> impedance_spectroscopy_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("impedance_spectroscopy"));
    std::fstream fout("computed_vs_exact_impedance_spectroscopy_data", std::fstream::out);
    for (double frequency = frequency_lower_limit; frequency < frequency_upper_limit; frequency *= ratio)
    {
        impedance_spectroscopy_database->put("frequency", frequency);
        std::complex<double> computed_impedance = cap::measure_impedance(device, impedance_spectroscopy_database);
        std::complex<double> exact_impedance    =
            ((type.compare("SeriesRC") == 0) ?
            series_resistance + 1.0 / std::complex<double>(0.0, capacitance * 2.0*pi*frequency) :
            series_resistance + parallel_resistance / std::complex<double>(1.0, capacitance * parallel_resistance * 2.0*pi*frequency));
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
          BOOST_CHECK_CLOSE(computed_impedance.real(), exact_impedance.real(), percent_tolerance);
          BOOST_CHECK_CLOSE(computed_impedance.imag(), exact_impedance.imag(), percent_tolerance);
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
    cap::bar(device, impedance_spectroscopy_database, fout);

}    
