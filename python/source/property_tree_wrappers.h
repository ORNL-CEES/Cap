#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/python/list.hpp>
#include <string>

namespace pycap {

using string = std::string;

double get_double(boost::property_tree::ptree const & ptree, string const & path);
string get_string(boost::property_tree::ptree const & ptree, string const & path);
int    get_int   (boost::property_tree::ptree const & ptree, string const & path);
bool   get_bool  (boost::property_tree::ptree const & ptree, string const & path);

double get_double_with_default_value(boost::property_tree::ptree const & ptree, string const & path, double const & default_value);
string get_string_with_default_value(boost::property_tree::ptree const & ptree, string const & path, string const & default_value);
int    get_int_with_default_value   (boost::property_tree::ptree const & ptree, string const & path, int    const & default_value);
bool   get_bool_with_default_value  (boost::property_tree::ptree const & ptree, string const & path, bool   const & default_value);

void put_double(boost::property_tree::ptree & ptree, string const & path, double const & value);
void put_string(boost::property_tree::ptree & ptree, string const & path, string const & value);
void put_int   (boost::property_tree::ptree & ptree, string const & path, int    const & value);
void put_bool  (boost::property_tree::ptree & ptree, string const & path, bool   const & value);

boost::python::list get_array_double(boost::property_tree::ptree const & ptree, string const & path);
boost::python::list get_array_string(boost::property_tree::ptree const & ptree, string const & path);
boost::python::list get_array_int   (boost::property_tree::ptree const & ptree, string const & path);
boost::python::list get_array_bool  (boost::property_tree::ptree const & ptree, string const & path);

void parse_xml (boost::property_tree::ptree & ptree, string const & filename);
void parse_json(boost::property_tree::ptree & ptree, string const & filename);

boost::property_tree::ptree get_child(boost::property_tree::ptree & ptree, string const & path);

}
