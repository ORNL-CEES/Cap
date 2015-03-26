#define BOOST_TEST_MODULE TestResistorCapacitor2   
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

BOOST_AUTO_TEST_CASE( test_resistor_capacitor )
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
    int    const ignore_cycles       = input_database->get<int   >("impedance_spectroscopy.ignore_cycles"  );
    int    const steps_per_cycle     = input_database->get<int   >("impedance_spectroscopy.steps_per_cycle");
    double const tolerance           = input_database->get<double>("impedance_spectroscopy.tolerance"      );
    double const initial_voltage     = 0.0;
    std::string const type = input_database->get<std::string>("device.type");


    double       time      = 0.0;
    double const time_step = 1.0 / frequency / steps_per_cycle;
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

    std::cout<<type<<"\n";
    std::cout<<"delta_t = "<<time_step<<"\n";
    double const tau = series_resistance*parallel_resistance/(series_resistance+parallel_resistance)*capacitance;
    std::cout<<"tau = "<<tau<<"\n";
    std::cout<<"delta_t / tau = "<<time_step/tau<<"\n";
    std::cout<<"expm1(- delta_t / tau) = "<<std::exp(-time_step/tau)<<"\n";

    double const angular_frequency = 2.0 * pi * frequency;

    double const gain  =
        (
            (type.compare("SeriesRC") == 0)
        ?
            angular_frequency * capacitance / std::sqrt(1.0 + std::pow(angular_frequency * series_resistance * capacitance, 2))
        :
            1.0 / (series_resistance + parallel_resistance) / (1.0 + std::pow(angular_frequency*series_resistance*parallel_resistance/(series_resistance+parallel_resistance)*capacitance,2))
                  * std::sqrt(
                      std::pow(1.0 + std::pow(angular_frequency * parallel_resistance * capacitance,2 ) * series_resistance / (series_resistance + parallel_resistance), 2) 
                    + std::pow(angular_frequency * std::pow(parallel_resistance, 2) / (series_resistance + parallel_resistance) * capacitance, 2)
                  )
        );

    double const phase = 
        (
            (type.compare("SeriesRC") == 0)
        ?
            std::atan(
                1.0 
                /
                (angular_frequency * series_resistance * capacitance)
            )
        :
            std::atan(
                angular_frequency * std::pow(parallel_resistance, 2) / (series_resistance + parallel_resistance) * capacitance
                /
                (1.0 + std::pow(angular_frequency * parallel_resistance * capacitance, 2) * series_resistance / (series_resistance + parallel_resistance))
            )
        );

    for (int n = 0; n < cycles*steps_per_cycle; ++n)
    {
          time += time_step;
          voltage = amplitude * std::sin(angular_frequency * time);
          device->evolve_one_time_step_constant_voltage(time_step, voltage);
          device->get_current(current);
          double const exact = amplitude * gain * std::sin(angular_frequency * time + phase);
          double const error = 100.0 * std::abs(current-exact)/(amplitude*gain);
          if (n >= ignore_cycles*steps_per_cycle)
              BOOST_CHECK_SMALL(error, tolerance);

          fout<<boost::format("  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  \n")
              % time
              % current
              % voltage
              % exact
              % error
              ;
    }

}    
