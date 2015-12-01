#define BOOST_TEST_MODULE ExactTransientSolution
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <cap/utils.h>
#include <tasmanian/TasmanianSparseGrid.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/math/distributions/beta.hpp>
#include <iostream>
#include <fstream>
#include <numeric>

namespace cap {

std::map<double,std::complex<double>>
run_code(std::shared_ptr<boost::property_tree::ptree const> database)
{
    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // measure its impedance
    std::shared_ptr<boost::property_tree::ptree> eis_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child("impedance_spectroscopy"));
    return measure_impedance(device, eis_database);
}



void do_stuff(std::shared_ptr<boost::property_tree::ptree> database)
{
    // beta distribution
    double const alpha  = database->get<double>("uq.alpha");
    double const beta   = database->get<double>("uq.beta" );
    int    const depth  = database->get<int   >("uq.depth");

    // read parameters from the database
    int    const params = database->get<int   >("uq.params");
    std::vector<double> mins(params);
    std::vector<double> maxs(params);
    std::vector<std::string> names(params);
    for (int p = 0; p < params; ++p)
    {
        boost::property_tree::ptree param_database =
            database->get_child("uq.param_"+std::to_string(p));
        std::string param_name = param_database.get<std::string>("name");
        std::vector<double> param_range =
            cap::to_vector<double>(param_database.get<std::string>("range"));
        std::string param_distribution_type =
            param_database.get<std::string>("distribution_type");
        BOOST_TEST(param_range.size() == 2);
        mins[p]  = param_range[0];
        maxs[p]  = param_range[1];
        names[p] = param_name;
    }
    // check that param exists in the database and that its nominal
    //     value is within the range that is specified
    for (int p = 0; p < params; ++p)
    {
        std::cout<<"param_"<<p<<"  "<<names[p]<<"  " 
            <<"[min max] = ["<<mins[p]<<" "<<maxs[p]<<"]\n";
        double val = std::nan("");
        // TODO: make it require no throw
        BOOST_CHECK_NO_THROW(
            val = database->get<double>(names[p]) );
        BOOST_TEST(val >= mins[p]);
        BOOST_TEST(val <= maxs[p]);
    }
    boost::math::beta_distribution<double> distribution(alpha, beta);

    TasGrid::TasmanianSparseGrid grid;
    grid.makeGlobalGrid(params, 2, depth, TasGrid::type_level, TasGrid::rule_gaussjacobi, 0, alpha, beta);
    grid.setDomainTransform(&(mins[0]), &(maxs[0]));
    std::cout<<grid.getNumPoints()<<"  "<<grid.getNumDimensions()+1<<"\n";
    std::cout<<"alpha = "<<alpha<<"\n";
    std::cout<<"beta  = "<<beta <<"\n";
    int    const   n_points = grid.getNumPoints        ();
    double const * weights  = grid.getQuadratureWeights();
    double const * points   = grid.getPoints           ();
    // check that the evalutation points are within the range
    for (int q = 0; q < n_points; ++q)
    {
        std::cout<<"  "<<weights[q]<<"  ";
        for (int p = 0; p < params; ++p)
            std::cout<<points[q*params+p]<<"  ";
        std::cout<<"\n";
        for (int p = 0; p < params; ++p)
        {
            BOOST_CHECK_GE(points[q*params+p], mins[p]);
            BOOST_CHECK_LE(points[q*params+p], maxs[p]);
        }
    }

    // run the code
    std::vector<std::map<double,std::complex<double>>> results;
    for (int q = 0; q < n_points; ++q)
    {
        for (int p = 0; p < params; ++p)
        {
            database->put(names[p], points[q*params+p]);
        }
        results.emplace_back(run_code(database));
        // print results
        std::cout<<"## "<<q<<" ##\n";
        for (auto x : results.back())
            std::cout<<x.first<<"  "<<x.second<<"\n";
    }

    std::map<double,std::complex<double>> expected_value;
    for (int q = 0; q < n_points; ++q)
    {
        for (auto x : results[q])
            expected_value[x.first] += weights[q] * x.second;
    }

    for (auto & x : expected_value)
        x.second /= std::accumulate(&(weights[0]), &(weights[n_points]), 0.0);

    std::cout<<"#######\n";
    for (auto x : expected_value)
        std::cout<<x.first<<"  "<<x.second<<"\n";

    // build an interpolant
    double frequency;
    std::complex<double> complex_impedance;
    std::vector<double> values(2*n_points);
std::cout<<"needed="<<grid.getNumNeeded()<<"\n";
std::cout<<"dimesnsions="<<grid.getNumDimensions()<<"\n";
std::cout<<"points="<<grid.getNumPoints()<<"\n";
std::cout<<"output="<<grid.getNumOutputs()<<"\n";
std::cout<<"loaded="<<grid.getNumLoaded()<<"\n";
    if (grid.getNumNeeded()*grid.getNumOutputs() != static_cast<int>(values.size()))
        throw std::runtime_error("check outputs number in the grid");
    for (auto x : expected_value)
    {
        std::tie(frequency, complex_impedance) = x;
        for (int q = 0; q < n_points; ++q)
        {
            values[2*q+0] = std::real(results[q][frequency]);
            values[2*q+1] = std::imag(results[q][frequency]);
        }
//        grid.loadNeededPoints(&(values[0]));
    }
double const dummy_frequency = expected_value.begin()->first;

    // compute the coefficients
    if (grid.getNumDimensions() != params)
        throw std::runtime_error("check dimensions number in the grid");
    std::vector<std::complex<double>> coefficients(n_points, 0.0);
    for (int q = 0; q < n_points; ++q)
    {
        double const * first = grid.getInterpolationWeights(&(points[q]));
        double const * last  = first+n_points;
        std::vector<double> interpolation_weights(first, last);
        for (auto h : interpolation_weights)
            std::cout<<h<<"  ";
        std::cout<<"\n";
        for (int p = 0; p < n_points; ++p)
            coefficients[p] += results[q][dummy_frequency] * interpolation_weights[p];
    }
    for (auto c : coefficients)
        std::cout<<c<<"\n";
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_uq_eis )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_uq_eis", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    cap::do_stuff(input_database);
}    
