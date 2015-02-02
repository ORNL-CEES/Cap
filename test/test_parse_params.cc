#define BOOST_TEST_MODULE MyFirstTest
#define BOOST_TEST_MAIN
#include <cap/utils.h>
#include <deal.II/base/types.h>
#include <deal.II/base/exceptions.h>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_parse_params )
{
   // create a database
   boost::property_tree::ptree params;

   // ensure material_ids and boundary_ids are parsed correctly
   params.put("some_material_id", 3);
   params.put("some_boundary_id", 1);
   BOOST_CHECK_EQUAL( params.get<dealii::types::material_id>("some_material_id"), dealii::types::material_id(3) );
   BOOST_CHECK_EQUAL( params.get<dealii::types::material_id>("some_boundary_id"), dealii::types::boundary_id(1) );

   // ensure property tree will throw
   params.put("is_a_string", "hello cruel world");
   BOOST_CHECK_THROW( params.get<double>("invalid_key"), boost::property_tree::ptree_bad_path );
   BOOST_CHECK_THROW( params.get<double>("is_a_string"), boost::property_tree::ptree_bad_data );

   // ensure to_string and to_vector are working properly
   params.put("ones",  cap::to_string(std::vector<double>     (10, 1.0)  ));
   params.put("yes",   cap::to_string(std::vector<std::string>(10, "yes")));
   params.put("empty", cap::to_string(std::vector<int>        ()         ));
   std::vector<double>      ones  = cap::to_vector<double>     (params.get<std::string>("ones") );
   std::vector<std::string> yes   = cap::to_vector<std::string>(params.get<std::string>("yes")  );
   std::vector<int>         empty = cap::to_vector<int>        (params.get<std::string>("empty"));
   BOOST_FOREACH(double const & val, ones) 
       BOOST_CHECK_EQUAL( val, 1.0 );
   BOOST_FOREACH(std::string const & val, yes) 
       BOOST_CHECK_EQUAL( val, "yes" );
   BOOST_CHECK( empty.empty() );
}
