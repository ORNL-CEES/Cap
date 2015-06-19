#define BOOST_TEST_MODULE TestResistorCapacitor2   
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>

#include <cap/resistor_capacitor.h>
#include <memory>

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
    BOOST_CHECK_EQUAL(std::acos(-1.0), boost::math::constants::pi<double>());

    device->reset_voltage(initial_voltage);
    double voltage;
    double current;
    std::fstream fout;
    fout.open("resistor_capacitor_data", std::fstream::out);

    std::cout<<type<<"\n";
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
          device->evolve_one_time_step_changing_voltage(time_step, voltage);
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



void foo(
    std::shared_ptr<cap::EnergyStorageDevice         >   step_device,
    std::shared_ptr<cap::EnergyStorageDevice         >   ramp_device,
    std::shared_ptr<boost::property_tree::ptree const>   database,
    std::ostream                                       & os
    )
{
    std::string const type              = database->get<std::string>("type"             );
    double      const series_resistance = database->get<double     >("series_resistance");
    double      const initial_voltage   = database->get<double     >("initial_voltage"  );
    double      const final_voltage     = database->get<double     >("final_voltage"    );
    double      const initial_time      = database->get<double     >("initial_time"     );
    double      const final_time        = database->get<double     >("final_time"       );
    int         const steps             = database->get<int        >("steps"            );
    double      const time_step         = (final_time - initial_time) / static_cast<double>(steps);

    double step_capacitor_voltage;
    double step_current;
    double step_voltage;
    double ramp_capacitor_voltage;
    double ramp_current;
    double ramp_voltage;
    for (double time = initial_time; time < final_time-time_step*std::numeric_limits<double>::epsilon(); time += time_step)
    {
        ramp_voltage = initial_voltage + (final_voltage - initial_voltage) / (final_time - initial_time) * (time+time_step - initial_time);
        step_voltage = final_voltage;
        ramp_device->evolve_one_time_step_changing_voltage(time_step, ramp_voltage);
        step_device->evolve_one_time_step_changing_voltage(time_step, step_voltage);
        if (type.compare("SeriesRC") == 0) {
            step_capacitor_voltage = (std::dynamic_pointer_cast<cap::SeriesRC  >(step_device))->U_C;
            ramp_capacitor_voltage = (std::dynamic_pointer_cast<cap::SeriesRC  >(ramp_device))->U_C;
        } else if (type.compare("ParallelRC") == 0) {
            step_capacitor_voltage = (std::dynamic_pointer_cast<cap::ParallelRC>(step_device))->U_C;
            ramp_capacitor_voltage = (std::dynamic_pointer_cast<cap::ParallelRC>(ramp_device))->U_C;
        } else {
            throw std::runtime_error("NEIN NEIN NEIN");
        }
        step_device->get_current(step_current);
        ramp_device->get_current(ramp_current);
        os<<boost::format("  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  \n")
            % (time+time_step)
            % step_current
            % step_voltage
            % step_capacitor_voltage
            % (series_resistance * step_current)
            % ramp_current
            % ramp_voltage
            % ramp_capacitor_voltage
            % (series_resistance * ramp_current)
            ;
    }
}



BOOST_AUTO_TEST_CASE( test_step_vs_ramp )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_resistor_capacitor", *input_database);

    // build two energy storage systems
    // one of them will be operated with steps
    // the one is operated as a ramp
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> step_device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));
    std::shared_ptr<cap::EnergyStorageDevice> ramp_device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    std::fstream fout;
    fout.open("step_vs_ramp_data", std::fstream::out);

    std::string const type = input_database->get<std::string>("device.type");
    double const series_resistance   = input_database->get<double>("device.series_resistance"              );
    double const parallel_resistance = input_database->get<double>("device.parallel_resistance"            );
    double const capacitance         = input_database->get<double>("device.capacitance"                    );
    double const initial_voltage     = input_database->get<double>("initial_voltage"                       );
    double const final_voltage       = input_database->get<double>("final_voltage"                         );
    double const time_constant       =
        (
            (type.compare("SeriesRC") == 0)
        ?
            series_resistance * capacitance
        :
            series_resistance * parallel_resistance / (series_resistance + parallel_resistance) * capacitance
        );

    step_device->reset_voltage(initial_voltage);
    ramp_device->reset_voltage(initial_voltage);
    double const initial_time       = input_database->get<double>("initial_time"                          , 0.0);
    double const final_time         = input_database->get<double>("final_time"                                 );
    int    const steps              = input_database->get<int   >("steps"                                      );
    int    const control_steps      = input_database->get<int   >("control_steps"                              );
    double const time_step          = (final_time - initial_time) / static_cast<double>(steps);

    std::cout<<"tau = "<<time_constant<<"\n";
    std::cout<<"delta_t = "<<final_time<<"\n";
    std::cout<<"delta_t/tau = "<<final_time/time_constant<<"\n";

    std::shared_ptr<boost::property_tree::ptree> foo_database =
        std::make_shared<boost::property_tree::ptree>();

    std::fstream fout2("step_vs_ramp_data2", std::fstream::out);
    foo_database->put("type"             , type             );
    foo_database->put("series_resistance", series_resistance);
    foo_database->put("steps"            , control_steps    );
    double voltage = initial_voltage;
    double step_current;
    double ramp_current;
    for (double time = initial_time; time < final_time; time += time_step)
    {
        foo_database->put("initial_time"     , time             );
        foo_database->put("final_time"       , time+time_step   );
        foo_database->put("initial_voltage"  , voltage          );
        voltage = 0.5 * (initial_voltage + final_voltage) + 0.5 * (final_voltage - initial_voltage) * boost::math::cos_pi((time+time_step - final_time) / (final_time - initial_time));
        std::cout<<time<<"  "<<voltage<<"\n";
        foo_database->put("final_voltage"    , voltage          );
        foo(step_device, ramp_device, foo_database, fout);
        step_device->get_current(step_current);
        ramp_device->get_current(ramp_current);
        fout2<<boost::format("  %22.15f  %22.15f  %22.15f  %22.15f  \n")
            % (time+time_step)
            % voltage
            % step_current
            % ramp_current
            ;
    }
    std::cout<<"step = "<<(series_resistance * step_current)<<"\n";
    std::cout<<"ramp = "<<(series_resistance * ramp_current)<<"\n";

    fout.close();

}
