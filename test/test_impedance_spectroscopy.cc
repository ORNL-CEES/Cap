#define BOOST_TEST_MODULE TestImpedanceSpectroscopy
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <iterator>
#include <algorithm>



namespace cap {

void impedance_spectroscopy(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double const frequency_upper_limit = database->get<double>("frequency_upper_limit");
    double const frequency_lower_limit = database->get<double>("frequency_lower_limit");
    int    const steps_per_decade      = database->get<int   >("steps_per_decade"     );
    double const pi                    = boost::math::constants::pi<double>();
    std::shared_ptr<boost::property_tree::ptree> tmp =
        std::make_shared<boost::property_tree::ptree>(*database);
    os<<"# impedance Z(f) = R + i X \n";
    os<<boost::format( "# %22s  %22s  %22s  %22s  %22s  \n")
        % "frequency_f_[Hz]"
        % "resistance_R_[ohm]"
        % "reactance_X_[ohm]"
        % "magnitude_|Z|_[ohm]"
        % "phase_arg(Z)_[degree]"
        ;
    double frequency;
    std::complex<double> impedance;
    for (double f= frequency_upper_limit; f>= frequency_lower_limit; f /= std::pow(10.0, 1.0/steps_per_decade))
    {
        tmp->put("frequency", f);
        for (auto const & data : measure_impedance(dev, tmp))
        {
            std::tie(frequency, impedance) = data;
            os<<boost::format( "  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  \n")
                % frequency
                % impedance.real()
                % impedance.imag()
                % std::abs(impedance)
                % (std::arg(impedance) * 180.0 / pi)
                ;
        }
    }
}

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_impedance_spectroscopy )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_impedance_spectroscopy", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

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

    fout.close();
}
