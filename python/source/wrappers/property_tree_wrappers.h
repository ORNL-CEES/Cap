#ifndef PROPERTY_TREE_WRAPPERS_H
#define PROPERTY_TREE_WRAPPERS_H

#include <boost/property_tree/ptree.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/object/pickle_support.hpp>
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
void parse_info(boost::property_tree::ptree & ptree, string const & filename);

boost::property_tree::ptree get_child(boost::property_tree::ptree & ptree, string const & path);
void put_child(boost::property_tree::ptree & ptree, string const & path, boost::property_tree::ptree const & child);

template <typename T>
struct serializable_class_pickle_support : boost::python::pickle_suite
{
    static
    boost::python::tuple
    getstate(T const & object)
    {
        std::ostringstream os;
        boost::archive::text_oarchive oa(os);
        oa << object;
        return boost::python::make_tuple(os.str());
    }

    static
    void
    setstate(T & object, boost::python::tuple state)
    {
        if (boost::python::len(state) != 1)
        {
          PyErr_SetObject(PyExc_ValueError,
                          ("expected 1-item tuple in call to __setstate__; got %s"
                           % state).ptr()
              );
          boost::python::throw_error_already_set();
        }

        boost::python::str pystr = boost::python::extract<boost::python::str>(state[0])();
        std::string cppstr = boost::python::extract<std::string>(pystr)();
        std::istringstream is(cppstr);
        boost::archive::text_iarchive ia(is);
        ia >> object;
    }
};

}

#endif
