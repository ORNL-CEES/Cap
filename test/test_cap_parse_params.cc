#include <cache/utils.h>
#include <deal.II/base/types.h>
#include <deal.II/base/exceptions.h>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>

int main(int argc, char *argv[])
{
   boost::property_tree::ptree params;

   // ensure material_ids and boundary_ids are parsed correctly
   params.put("some_material_id", 3);
   params.put("some_boundary_id", 1);
   Assert(params.get<dealii::types::material_id>("some_material_id") == dealii::types::material_id(3),
       dealii::StandardExceptions::ExcInternalError());
   Assert(params.get<dealii::types::boundary_id>("some_boundary_id") == dealii::types::boundary_id(1),
       dealii::StandardExceptions::ExcInternalError());

   // ensure to_string and to_vector are working properly
   params.put("ones", cache::to_string(std::vector<double>(10, 1.0)));
   params.put("yes",  cache::to_string(std::vector<std::string>(10, "yes")));
   std::vector<double> ones = cache::to_vector<double>(params.get<std::string>("ones"));
   std::vector<std::string> yes = cache::to_vector<std::string>(params.get<std::string>("yes"));
   BOOST_FOREACH(double const & val, ones) 
       Assert(val == 1.0, dealii::StandardExceptions::ExcInternalError());
   BOOST_FOREACH(std::string const & val, yes) 
       Assert(val == "yes", dealii::StandardExceptions::ExcInternalError());

   return 0;
}
