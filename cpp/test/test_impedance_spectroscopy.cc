/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#define BOOST_TEST_MODULE TestImpedanceSpectroscopy
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>



namespace cap {

void run_eis(std::shared_ptr<cap::EnergyStorageDevice> device,
    std::shared_ptr<boost::property_tree::ptree const> eis_database,
    std::ostream & os = std::cout)
{
    // run the impedance spectroscopy experiment
    std::map<double, std::complex<double>> eis_data =
        impedance_spectroscopy(device, eis_database);

    // report the results
    double const pi = boost::math::constants::pi<double>();
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
    for (auto const & data : eis_data)
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

} // end namespace cap



BOOST_AUTO_TEST_CASE( test_impedance_spectroscopy )
{
    // parse input file
    boost::property_tree::ptree input_database;
    boost::property_tree::xml_parser::read_xml("input_impedance_spectroscopy", input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    // build an energy storage device
    boost::property_tree::ptree device_database =
        input_database.get_child("device");
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(boost::mpi::communicator(), device_database);

    // measure its impedance
    std::fstream fout;
    fout.open("impedance_spectroscopy_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> impedance_spectroscopy_database =
        std::make_shared<boost::property_tree::ptree>(input_database.get_child("impedance_spectroscopy"));
    cap::run_eis(device, impedance_spectroscopy_database, fout);

    fout.close();
}
