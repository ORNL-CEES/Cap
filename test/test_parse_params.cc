#define BOOST_TEST_MODULE MyFirstTest
#define BOOST_TEST_MAIN
#include <cap/utils.h>
#ifdef WITH_DEAL_II
  #include <deal.II/base/types.h>
  #include <deal.II/base/exceptions.h>
#endif
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_parse_params )
{
   // create a database
   boost::property_tree::ptree params;

#ifdef WITH_DEAL_II
   // ensure material_ids and boundary_ids are parsed correctly
   params.put("some_material_id", 3);
   params.put("some_boundary_id", 1);
   BOOST_CHECK_EQUAL( params.get<dealii::types::material_id>("some_material_id"), dealii::types::material_id(3) );
   BOOST_CHECK_EQUAL( params.get<dealii::types::material_id>("some_boundary_id"), dealii::types::boundary_id(1) );

   params.put("material_ids", "0,1,2,3");
   std::vector<dealii::types::material_id> material_ids = cap::to_vector<dealii::types::material_id>(params.get<std::string>(("material_ids")));
   BOOST_CHECK_EQUAL( material_ids[0], dealii::types::material_id(0) );
   BOOST_CHECK_EQUAL( material_ids[1], dealii::types::material_id(1) );
   BOOST_CHECK_EQUAL( material_ids[2], dealii::types::material_id(2) );
   BOOST_CHECK_EQUAL( material_ids[3], dealii::types::material_id(3) );
#endif

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
   std::vector<bool> vec_in = { true, false, false };
   params.put("boule", cap::to_string(vec_in));
   std::vector<bool> vec_out = cap::to_vector<bool>(params.get<std::string>("boule"));
   BOOST_CHECK_EQUAL_COLLECTIONS(vec_in.begin(), vec_in.end(), vec_out.begin(), vec_out.end());

   params.put("always_true", "    true,true, true  ");
   std::vector<bool> always_true = cap::to_vector<bool>(params.get<std::string>("always_true"));
   for (bool const is_true : always_true)
       BOOST_TEST( is_true );
   
   // space allowed
   std::vector<double> bugfixme = { 1.234, 5.678 };
   params.put("bugfixme", "1.234 ,    5.678 ");
   std::vector<double> fixed =
       cap::to_vector<double>(params.get<std::string>("bugfixme"));
   BOOST_TEST(bugfixme == fixed);

   // ensure entry is overwritten
   std::string const default_value("default");
   std::string const path("some.path");
   params.put(path, default_value);
   BOOST_CHECK_EQUAL(default_value.compare(params.get<std::string>(path)), 0);
   std::string const new_value("new");
   params.put(path, new_value);
   BOOST_CHECK_EQUAL(new_value.compare(params.get<std::string>(path)), 0);
   BOOST_CHECK_THROW(params.get<double>(path), boost::property_tree::ptree_bad_data);
   double const pi = 3.14159265359;
   params.put(path, pi);
   BOOST_CHECK_EQUAL(params.get<double>(path), pi);
   std::string const pi_to_string = std::to_string(pi);
   std::size_t const precision = 6;
   BOOST_CHECK_EQUAL(pi_to_string.compare(0, precision+1, params.get<std::string>(path), 0, precision+1), 0);
    
   
}
