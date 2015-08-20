#define BOOST_TEST_MODULE ExactTransientSolution
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/utils.h>
#include <tasmanian/TasmanianSparseGrid.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/math/distributions/beta.hpp>
#include <iostream>
#include <fstream>

namespace cap {

double run_code(std::shared_ptr<boost::property_tree::ptree> database)
{
    return std::nan("");
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
    grid.makeGlobalGrid(params, 0, depth, TasGrid::type_level, TasGrid::rule_gaussjacobi, 0, alpha, beta);
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
        // run the code
    }
    std::vector<double> results;
    for (int q = 0; q < n_points; ++q)
    {
        for (int p = 0; p < params; ++p)
        {
            database->put(names[p], points[q*params+p]);
        }
        results.emplace_back(run_code(database));
   }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_uq_eis )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_uq_eis", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    cap::do_stuff(input_database);
}    
