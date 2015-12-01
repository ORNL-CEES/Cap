#include <cap/version.h>
#include <pycap/property_tree_wrappers.h>
#include <pycap/energy_storage_device_wrappers.h>
#include <boost/python.hpp>
#include <string>
#include <memory>
#include <map>


BOOST_PYTHON_MODULE(_pycap)
{
    boost::python::class_<pycap::ElectrochemicalImpedanceSpectroscopyData, std::shared_ptr<pycap::ElectrochemicalImpedanceSpectroscopyData>>("ElectrochemicalImpedanceSpectroscopyData")
        .def("impedance_spectroscopy", &pycap::ElectrochemicalImpedanceSpectroscopyData::impedance_spectroscopy)
        .def("measure_impedance", &pycap::ElectrochemicalImpedanceSpectroscopyData::measure_impedance)
        .def("get_frequency", &pycap::ElectrochemicalImpedanceSpectroscopyData::get_frequency)
        .def("get_complex_impedance", &pycap::ElectrochemicalImpedanceSpectroscopyData::get_complex_impedance)
        .def("clear", &pycap::ElectrochemicalImpedanceSpectroscopyData::clear)
        ;
    boost::python::class_<pycap::EnergyStorageDeviceWrap, std::shared_ptr<pycap::EnergyStorageDeviceWrap>, boost::noncopyable>("EnergyStorageDevice", "Wrappers for Cap.EnergyStorageDevice", boost::python::no_init)
        .def("__init__", boost::python::make_constructor(&pycap::build_energy_storage_device) )
        .def("get_voltage", (&pycap::get_voltage), boost::python::args("self") )
        .def("get_current", (&pycap::get_current), boost::python::args("self") )
        .def("evolve_one_time_step_constant_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_current), boost::python::args("self", "time_step", "current") )
        .def("evolve_one_time_step_constant_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_voltage), boost::python::args("self", "time_step", "voltage") )
        .def("evolve_one_time_step_constant_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_power  ), boost::python::args("self", "time_step", "load"   ) )
        .def("evolve_one_time_step_constant_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_load   ), boost::python::args("self", "time_step", "power"  ) )
        .def("evolve_one_time_step_changing_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_current), boost::python::args("self", "time_step", "current") )
        .def("evolve_one_time_step_changing_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_voltage), boost::python::args("self", "time_step", "voltage") )
        .def("evolve_one_time_step_changing_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_power  ), boost::python::args("self", "time_step", "power"  ) )
        .def("evolve_one_time_step_changing_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_load   ), boost::python::args("self", "time_step", "load"   ) )
        .def("compute_equivalent_circuit", &pycap::compute_equivalent_circuit, "Return the PropertyTree to build an equivalent circuit model", boost::python::args("ptree") )
        .staticmethod("compute_equivalent_circuit")
        ;
    boost::python::register_ptr_to_python<std::shared_ptr<cap::EnergyStorageDevice>>();

    boost::python::scope().attr("__version__"        ) = cap::version()        ;
    boost::python::scope().attr("__git_branch__"     ) = cap::git_branch()     ;
    boost::python::scope().attr("__git_commit_hash__") = cap::git_commit_hash();

    boost::python::class_<boost::property_tree::ptree,std::shared_ptr<boost::property_tree::ptree>>("PropertyTree", "Wrappers for Boost.PropertyTree")
        .def("get_double"                   , &pycap::get_double                   , "Get the double at the given path."                                , boost::python::args("self", "path") )
        .def("get_string"                   , &pycap::get_string                   , "Get the string at the given path."                                , boost::python::args("self", "path") )
        .def("get_int"                      , &pycap::get_int                      , "Get the integer at the given path."                               , boost::python::args("self", "path") )
        .def("get_bool"                     , &pycap::get_bool                     , "Get the boolean at the given path."                               , boost::python::args("self", "path") )
        .def("get_double_with_default_value", &pycap::get_double_with_default_value, "Get the double at the given path or return default_value."        , boost::python::args("self", "path", "default_value") )
        .def("get_string_with_default_value", &pycap::get_string_with_default_value, "Get the string at the given path or return default_value."        , boost::python::args("self", "path", "default_value") )
        .def("get_int_with_default_value"   , &pycap::get_int_with_default_value   , "Get the integer at the given path or return default_value."       , boost::python::args("self", "path", "default_value") )
        .def("get_bool_with_default_value"  , &pycap::get_bool_with_default_value  , "Get the boolean at the given path or return default_value."       , boost::python::args("self", "path", "default_value") )
        .def("put_double"                   , &pycap::put_double                   , "Set the node at the given path to the supplied value."            , boost::python::args("self", "path", "value") )
        .def("put_string"                   , &pycap::put_string                   , "Set the node at the given path to the supplied value."            , boost::python::args("self", "path", "value") )
        .def("put_int"                      , &pycap::put_int                      , "Set the node at the given path to the supplied value."            , boost::python::args("self", "path", "value") )
        .def("put_bool"                     , &pycap::put_bool                     , "Set the node at the given path to the supplied value."            , boost::python::args("self", "path", "value") )
        .def("get_array_double"             , &pycap::get_array_double             , "Get comma separated array of double."                             , boost::python::args("self", "path") )
        .def("get_array_string"             , &pycap::get_array_string             , "Get comma separated array of string."                             , boost::python::args("self", "path") )
        .def("get_array_int"                , &pycap::get_array_int                , "Get comma separated array of integer."                            , boost::python::args("self", "path") )
        .def("get_array_bool"               , &pycap::get_array_bool               , "Get comma separated array of boolean."                            , boost::python::args("self", "path") )
        .def("parse_xml"                    , &pycap::parse_xml                    , "Read the input file at XML format and populate the PropertyTree." , boost::python::args("self", "filename") )
        .def("parse_json"                   , &pycap::parse_json                   , "Read the input file at JSON format and populate the PropertyTree.", boost::python::args("self", "filename") )
        .def("get_child"                    , &pycap::get_child                    , "Get the child at the given path, or throw ptree_bad_path."        , boost::python::args("self", "path") )
        ;
}
