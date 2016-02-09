/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#define BOOST_TEST_MODULE ButlerVolmerKinetics
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <cmath>
#include <iostream>
#include <fstream>

BOOST_AUTO_TEST_CASE( test_butler_volmer_kinetic )
{
    double const faraday_constant                     =   9.64853365e4;
    double const gas_constant                         =   8.3144621;
    double const boltzmann_constant                   =   8.6173324e-5;
    double const temperature                          = 300.0;
    double const anodic_charge_transfer_coefficient   =   0.5;
    double const cathodic_charge_transfer_coefficient =   0.5;
    double const exchange_current_density             =   1.0;
    BOOST_TEST(
        faraday_constant/(gas_constant*temperature) ==
        1.0/(boltzmann_constant*temperature),
        boost::test_tools::tolerance(1.0e-6));
    std::cout<<boost::format("  %15.7e  %15.7e  \n")
        % static_cast<double>(faraday_constant/(gas_constant*temperature))
        % static_cast<double>(1.0/(boltzmann_constant*temperature))
        ;



    std::fstream fout;
    fout.open("butler_volmer_kinetics_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_butler_volmer_kinetics", *database);

    double const overpotential_lower_limit = database->get<double>("overpotential_lower_limit");
    double const overpotential_upper_limit = database->get<double>("overpotential_upper_limit");
    int    const n_points                  = database->get<int   >("n_points"                 );
    double const step = (overpotential_upper_limit - overpotential_lower_limit) / (n_points + 1);
    for (double overpotential = overpotential_lower_limit;
        overpotential <= overpotential_upper_limit;
        overpotential += step)
    {
        fout<<boost::format("  %10.7f  %15.7e  %15.7e  \n")
            % overpotential
            % static_cast<double>(
                exchange_current_density * 
                (std::exp(anodic_charge_transfer_coefficient * overpotential / (boltzmann_constant*temperature))
                - std::exp(- cathodic_charge_transfer_coefficient * overpotential / (boltzmann_constant*temperature)))
              )
            % static_cast<double>(
                exchange_current_density *
                (anodic_charge_transfer_coefficient + cathodic_charge_transfer_coefficient)
                    * overpotential / (boltzmann_constant*temperature)
              )
            ;

    }

    fout.close();
}
