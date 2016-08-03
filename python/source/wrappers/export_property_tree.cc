#include <pycap/property_tree_wrappers.h>
#include <boost/python/args.hpp>
#include <boost/python/class.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/exception_translator.hpp>

namespace pycap
{

char const * property_tree_docstring =
  "Wrappers for Boost.PropertyTree                                        \n"
  "                                                                       \n"
  "Examples                                                               \n"
  "--------                                                               \n"
  ">>> ptree = PropertyTree()                                             \n"
  ">>> ptree.put_double('pi', 3.14)                                       \n"
  ">>> ptree.get_double('pi')                                             \n"
  "3.14                                                                   \n"
  ">>> ptree.get_double_with_default('sqrt2', 1,41)                       \n"
  "1.41                                                                   \n"
  "                                                                       \n"
  "Raises                                                                 \n"
  "------                                                                 \n"
  "KeyError: No such node (<path>)                                        \n"
  "    Error indicating that specified <path> does not exist.             \n"
  "TypeError: conversion of data to type \"<type>\" failed                \n"
  "    Error indicating that translation from or to <type> has failed.    \n"
  "                                                                       \n"
  ;

std::string get_docstring(std::string const & type)
{
  return "Get the " + type + " at the given path.";
}

std::string get_with_default_value_docstring(std::string const & type)
{
  return "Get the " + type + " at the given path or return default_value.";
}

std::string put_docstring()
{
  return "Set the node at the given path to the supplied value.";
}

std::string get_array_docstring(std::string const & type)
{
  return "Get comma separated array of " + type + ".";
}

std::string parse_docstring(std::string const & format)
{
  return "Read the input file at " + format + " format and populate the PropertyTree.";
}

void export_property_tree()
{
  boost::python::class_<boost::property_tree::ptree>(
    "PropertyTree",
    property_tree_docstring)
    .def("get_double", &pycap::get_double, get_docstring("double").c_str(), boost::python::args("self", "path") )
    .def("get_string", &pycap::get_string, get_docstring("string").c_str(), boost::python::args("self", "path") )
    .def("get_int"   , &pycap::get_int   , get_docstring("int"   ).c_str(), boost::python::args("self", "path") )
    .def("get_bool"  , &pycap::get_bool  , get_docstring("bool"  ).c_str(), boost::python::args("self", "path") )
    .def("get_double_with_default_value", &pycap::get_double_with_default_value, get_with_default_value_docstring("double").c_str(), boost::python::args("self", "path", "default_value") )
    .def("get_string_with_default_value", &pycap::get_string_with_default_value, get_with_default_value_docstring("string").c_str(), boost::python::args("self", "path", "default_value") )
    .def("get_int_with_default_value"   , &pycap::get_int_with_default_value   , get_with_default_value_docstring("int"   ).c_str(), boost::python::args("self", "path", "default_value") )
    .def("get_bool_with_default_value"  , &pycap::get_bool_with_default_value  , get_with_default_value_docstring("bool"  ).c_str(), boost::python::args("self", "path", "default_value") )
    .def("put_double", &pycap::put_double, put_docstring().c_str(), boost::python::args("self", "path", "value") )
    .def("put_string", &pycap::put_string, put_docstring().c_str(), boost::python::args("self", "path", "value") )
    .def("put_int"   , &pycap::put_int   , put_docstring().c_str(), boost::python::args("self", "path", "value") )
    .def("put_bool"  , &pycap::put_bool  , put_docstring().c_str(), boost::python::args("self", "path", "value") )
    .def("get_array_double", &pycap::get_array_double, get_array_docstring("double").c_str(), boost::python::args("self", "path") )
    .def("get_array_string", &pycap::get_array_string, get_array_docstring("string").c_str(), boost::python::args("self", "path") )
    .def("get_array_int"   , &pycap::get_array_int   , get_array_docstring("int"   ).c_str(), boost::python::args("self", "path") )
    .def("get_array_bool"  , &pycap::get_array_bool  , get_array_docstring("bool"  ).c_str(), boost::python::args("self", "path") )
    .def("parse_xml" , &pycap::parse_xml , parse_docstring("XML" ).c_str(), boost::python::args("self", "filename") )
    .def("parse_json", &pycap::parse_json, parse_docstring("JSON").c_str(), boost::python::args("self", "filename") )
    .def("parse_info", &pycap::parse_info, parse_docstring("INFO").c_str(), boost::python::args("self", "filename") )
    .def("get_child", &pycap::get_child, "Get the child at the given path, or throw ptree_bad_path.", boost::python::args("self", "path") )
    .def("put_child", &pycap::put_child, "Put the child at the given path, create any missing parents, replace if it already exists.", boost::python::args("self", "path", "PropertyTree") )
    .def_pickle(pycap::serializable_class_pickle_support<boost::property_tree::ptree>())
        ;
  boost::python::register_exception_translator<boost::property_tree::ptree_error>(&pycap::translate);
  boost::python::register_ptr_to_python<std::shared_ptr<boost::property_tree::ptree>>();
}

} // end namespace pycap
