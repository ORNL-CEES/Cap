/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <pycap/property_tree_wrappers.h>
#include <cap/utils.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <vector>

namespace pycap {

double get_double(boost::property_tree::ptree const & ptree, string const & path)
{
    return ptree.get<double>(path);
}

string get_string(boost::property_tree::ptree const & ptree, string const & path)
{
    return ptree.get<string>(path);
}

int    get_int   (boost::property_tree::ptree const & ptree, string const & path)
{
    return ptree.get<int>(path);
}

bool   get_bool  (boost::property_tree::ptree const & ptree, string const & path)
{
    return ptree.get<bool>(path);
}

double get_double_with_default_value(boost::property_tree::ptree const & ptree, string const & path, double const & default_value)
{
    return ptree.get<double>(path, default_value);
}

string get_string_with_default_value(boost::property_tree::ptree const & ptree, string const & path, string const & default_value)
{
    return ptree.get<string>(path, default_value);
}

int    get_int_with_default_value   (boost::property_tree::ptree const & ptree, string const & path, int    const & default_value)
{
    return ptree.get<int>(path, default_value);
}

bool   get_bool_with_default_value  (boost::property_tree::ptree const & ptree, string const & path, bool   const & default_value)
{
    return ptree.get<bool>(path, default_value);
}

void put_double(boost::property_tree::ptree & ptree, string const & path, double const & value)
{
    ptree.put(path, value);
}

void put_string(boost::property_tree::ptree & ptree, string const & path, string const & value)
{
    ptree.put(path, value);
}

void put_int   (boost::property_tree::ptree & ptree, string const & path, int    const & value)
{
    ptree.put(path, value);
}

void put_bool  (boost::property_tree::ptree & ptree, string const & path, bool   const & value)
{
    ptree.put(path, value);
}

template <typename T>
boost::python::list get_array(boost::property_tree::ptree const & ptree, string const & path)
{
    std::vector<T> vector = 
        cap::to_vector<T>(ptree.get<string>(path));
    boost::python::list list;
    for (auto const & x : vector)
        list.append(x);
    return list;
}

boost::python::list get_array_double(boost::property_tree::ptree const & ptree, string const & path)
{
    return get_array<double>(ptree, path);
}

boost::python::list get_array_string(boost::property_tree::ptree const & ptree, string const & path)
{
    return get_array<string>(ptree, path);
}

boost::python::list get_array_int   (boost::property_tree::ptree const & ptree, string const & path)
{
    return get_array<int   >(ptree, path);
}

// the get_array function won't work with bool
// I believe it has to do with vector<bool>::operator[] returning a proxy object
boost::python::list get_array_bool  (boost::property_tree::ptree const & ptree, string const & path)
{
    std::vector<bool> vector = 
        cap::to_vector<bool>(ptree.get<string>(path));
    boost::python::list list;
    for (bool const x : vector)
        list.append(x);
    return list;
}

void parse_xml (boost::property_tree::ptree & ptree, string const & filename)
{
    boost::property_tree::xml_parser::read_xml(filename, ptree,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);
}

void parse_json(boost::property_tree::ptree & ptree, string const & filename)
{
    boost::property_tree::json_parser::read_json(filename, ptree);
}

void parse_info(boost::property_tree::ptree & ptree, string const & filename)
{
    boost::property_tree::info_parser::read_info(filename, ptree);
}

boost::property_tree::ptree get_child(boost::property_tree::ptree & ptree, string const & path)
{
    return ptree.get_child(path);
}

void put_child(boost::property_tree::ptree & ptree, string const & path, boost::property_tree::ptree const & child)
{
    ptree.put_child(path, child);
}

} // end namespace pycap
