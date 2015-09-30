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
    boost::python::class_<pycap::EnergyStorageDeviceWrap, std::shared_ptr<pycap::EnergyStorageDeviceWrap>, boost::noncopyable>("EnergyStorageDevice", boost::python::no_init)
        .def("__init__", boost::python::make_constructor(&pycap::build_energy_storage_device) )
        .def("get_voltage", (&pycap::get_voltage) )
        .def("get_current", (&pycap::get_current) )
        .def("evolve_one_time_step_constant_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_current) )
        .def("evolve_one_time_step_constant_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_voltage) )
        .def("evolve_one_time_step_constant_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_power  ) )
        .def("evolve_one_time_step_constant_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_load   ) )
        .def("evolve_one_time_step_changing_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_current) )
        .def("evolve_one_time_step_changing_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_voltage) )
        .def("evolve_one_time_step_changing_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_power  ) )
        .def("evolve_one_time_step_changing_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_load   ) )
        .def("compute_equivalent_circuit", &pycap::compute_equivalent_circuit)
        .staticmethod("compute_equivalent_circuit")
        ;
    boost::python::register_ptr_to_python<std::shared_ptr<cap::EnergyStorageDevice>>();

    boost::python::scope().attr("__version__"        ) = cap::version()        ;
    boost::python::scope().attr("__git_branch__"     ) = cap::git_branch()     ;
    boost::python::scope().attr("__git_commit_hash__") = cap::git_commit_hash();

    boost::python::class_<boost::property_tree::ptree, std::shared_ptr<boost::property_tree::ptree>>("PropertyTree")
        .def("get_double"                   , &pycap::get_double                   )
        .def("get_string"                   , &pycap::get_string                   )
        .def("get_int"                      , &pycap::get_int                      )
        .def("get_bool"                     , &pycap::get_bool                     )
        .def("get_double_with_default_value", &pycap::get_double_with_default_value)
        .def("get_string_with_default_value", &pycap::get_string_with_default_value)
        .def("get_int_with_default_value"   , &pycap::get_int_with_default_value   )
        .def("get_bool_with_default_value"  , &pycap::get_bool_with_default_value  )
        .def("put_double"                   , &pycap::put_double                   )
        .def("put_string"                   , &pycap::put_string                   )
        .def("put_int"                      , &pycap::put_int                      )
        .def("put_bool"                     , &pycap::put_bool                     )
        .def("get_array_double"             , &pycap::get_array_double             )
        .def("get_array_string"             , &pycap::get_array_string             )
        .def("get_array_int"                , &pycap::get_array_int                )
        .def("get_array_bool"               , &pycap::get_array_bool               )
        .def("parse_xml"                    , &pycap::parse_xml                    )
        .def("parse_json"                   , &pycap::parse_json                   )
        .def("get_child"                    , &pycap::get_child                    )
        ;
}
