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

BOOST_AUTO_TEST_CASE( test_impedance_spectroscopy )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_resistor_capacitor", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));


    double const series_resistance   = input_database->get<double>("device.series_resistance"              );
    double const parallel_resistance = input_database->get<double>("device.parallel_resistance"            );
    double const capacitance         = input_database->get<double>("device.capacitance"                    );
    double const frequency           = input_database->get<double>("impedance_spectroscopy.frequency"      );
    double const amplitude           = input_database->get<double>("impedance_spectroscopy.amplitude"      );
    int    const cycles              = input_database->get<int   >("impedance_spectroscopy.cycles"         );
    int    const steps_per_cycle     = input_database->get<int   >("impedance_spectroscopy.steps_per_cycle");
    double const initial_voltage     = 0.0;
    std::string const type = input_database->get<std::string>("device.type");


    double       time      = 0.0;
    double const time_step = 1.0 / frequency / steps_per_cycle;
    double const phase     = std::asin(initial_voltage / amplitude);
    double const pi        = std::acos(-1.0);
//    std::cout<<boost::format("  %20.15f  %20.15f  \n")
//        % std::acos(-1.0)
//        % boost::math::constants::pi<double>();
    BOOST_CHECK_EQUAL(std::acos(-1.0), boost::math::constants::pi<double>());

    device->reset_voltage(initial_voltage);
    double voltage;
    double current;
    std::fstream fout;
    fout.open("resistor_capacitor_data", std::fstream::out);

    std::cout<<"delta_t = "<<time_step<<"\n";
    std::cout<<"tau = "<<series_resistance*capacitance<<"\n";
    std::cout<<"delta_t / tau = "<<time_step/(series_resistance*capacitance)<<"\n";
    std::cout<<"exp(- delta_t / tau) = "<<std::exp(-time_step/(series_resistance*capacitance))<<"\n";

    double const harmonic_resistance =
        series_resistance * parallel_resistance / (series_resistance + parallel_resistance);
    double const angular_frequency =
        2.0 * pi * frequency;

//    for (int n = 0; n < cycles*steps_per_cycle; ++n)
    for (; time <= cycles / frequency; )
    {
          time += time_step;
          voltage = amplitude*std::sin(2.0*pi*frequency*time+phase);
          device->evolve_one_time_step_constant_voltage(time_step, voltage);
          device->get_current(current);
          double exact =
              ((type.compare("SeriesRC") == 0)  ?
              amplitude * 2.0*pi*frequency*capacitance / std::sqrt(1.0 + std::pow(2.0*pi*frequency*series_resistance*capacitance,2))
                  * std::sin(2.0*pi*frequency*time+phase+std::atan(1.0/(2.0*pi*frequency*series_resistance*capacitance)))
              :
              amplitude / (series_resistance+parallel_resistance) / (1.0 + std::pow(2.0*pi*frequency*series_resistance*parallel_resistance/(series_resistance+parallel_resistance)*capacitance,2))
                  * std::sqrt(std::pow(1.0+std::pow(2.0*pi*frequency*parallel_resistance*capacitance,2)*series_resistance/(series_resistance+parallel_resistance),2) + std::pow(2.0*pi*frequency*std::pow(parallel_resistance,2)/(series_resistance+parallel_resistance)*capacitance,2))
                  * std::sin(2.0*pi*frequency*time+phase+std::atan(2.0*pi*frequency*std::pow(parallel_resistance,2)/(series_resistance+parallel_resistance*capacitance)/(1.0+std::pow(2.0*pi*frequency*parallel_resistance*capacitance,2)*series_resistance/(series_resistance+parallel_resistance))))
              );

          fout<<boost::format("  %22.15e  %22.15e  %22.15e  %22.15e  \n")
              % time
              % current
              % voltage
              % exact
              ;
    }

}    
