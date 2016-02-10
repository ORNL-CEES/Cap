/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/version.h>
#include <pycap/property_tree_wrappers.h>
#include <pycap/energy_storage_device_wrappers.h>
#include <boost/python.hpp>
#include <string>
#include <memory>
#include <map>

#include <mpi4py/mpi4py.h>

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & py_ptree,
                            boost::python::object & py_comm)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(py_ptree);
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
    if (comm_p == NULL) boost::python::throw_error_already_set();
    boost::mpi::communicator comm(*comm_p, boost::mpi::comm_attach);
    return cap::EnergyStorageDevice::build(comm, ptree);
}

BOOST_PYTHON_MODULE(_pycap)
{
    if (import_mpi4py() < 0) return;

    boost::python::scope().attr("__version__"        ) = cap::version()        ;
    boost::python::scope().attr("__git_branch__"     ) = cap::git_branch()     ;
    boost::python::scope().attr("__git_commit_hash__") = cap::git_commit_hash();
    boost::python::scope().attr("__git_remote_url__" ) = cap::git_remote_url() ;
    boost::python::scope().attr("__doc__") =
        "\n"
        "PyCap\n"
        "=====\n"
        "Provides\n"
        "  1. Energy storage device models\n"
        "      Cap comes with a finite element model for supercapacitors as well as\n"
        "      equivalents circuits models.\n"
        "  2. Electrochemical measurement techniques\n"
        "      Cap can run a wide range of techniques on the available energy storage devices,\n"
        "      from the simple recording of discharge curves to the more complex calculation \n"
        "      of impedance spectra.\n"
        "\n"
        "How to use the documentation\n"
        "----------------------------\n"
        "Documentation is available in two forms: docstrings provided\n"
        "with the code, and a loose standing reference guide, available from\n"
        "`Cap's online documentation <https://cap.readthedocs.org>`_.\n"
        "\n"
        "Use the built-in ``help`` function to view a function's docstring::\n"
        "  >>> import pycap\n"
        "  >>> help(pycap.PropertyTree)\n"
        "  ...\n"
        "\n"
        "Utilities \n"
        "----------\n"
        "PropertyTree\n"
        "    Wrappers for Boost.PropertyTree\n"
        "    See documentation\n"
        "\n"
        "Available energy storage devices\n"
        "--------------------------------\n"
        "EnergyStorageDevice\n"
        "    Wrappers for Cap.EnergyStorageDevice\n"
        "    See documentation\n"
        "\n"
        "Available electrochemical techniques\n"
        "------------------------------------\n"
        "Discharge\n"
        "    Performs the discharge in one of four different control modes\n"
        "Charge\n"
        "    Performs the charge with an optional constant voltage step at the end\n"
        "CyclicChargeDischarge\n"
        "    Records charge and discharge curves through a number of cycles\n"
        "CyclicVoltammetry\n"
        "    Applies cyclic linear voltage ramps\n"
        "measure_impedance_spectrum\n"
        "    Measures the complex impedance over a range of frequencies\n"
        "measure_performance\n"
        "    Complete a series of discharges at various rate to produce a ragone plot\n"
        ;
    boost::python::docstring_options doc_options;
    doc_options.enable_user_defined();
    doc_options.enable_py_signatures();
    doc_options.disable_cpp_signatures();

    boost::python::class_<
        pycap::EnergyStorageDeviceWrap,
        std::shared_ptr<pycap::EnergyStorageDeviceWrap>,
        boost::noncopyable >(
        "EnergyStorageDevice",
        "Wrappers for Cap.EnergyStorageDevice                                \n"
        "                                                                    \n"
        "Examples                                                            \n"
        "--------                                                            \n"
        "                                                                    \n"
        ">>> ptree = PropertyTree()                                          \n"
        ">>> ptree.parse_info('device.info')                                 \n"
        ">>> device = EnergyStorageDevice(ptree)                             \n"
        ">>> delta_t, U = 1.0, 2.1 # units of seconds and volts              \n"
        ">>> device.evolve_one_time_step_constant_voltage(delta_t,U)         \n"
        ">>> I = device.get_current() # amperes                              \n"
        "                                                                    \n"
        ,
        boost::python::no_init )
        .def("__init__",
        boost::python::make_constructor(build_energy_storage_device,
//        boost::python::make_constructor(&pycap::build_energy_storage_device,
        boost::python::default_call_policies(),
        (boost::python::arg("ptree"), boost::python::arg("comm")="MPI.COMM_WORLD")),
        "                                                                    \n"
        "Parameters                                                          \n"
        "----------                                                          \n"
        "ptree : pycap.PropertyTree                                          \n"
        "    The appropriate property tree to create a device                \n"
        "    from the factory.                                               \n"
        "comm  : mpi4py.MPI.Comm                                             \n"
        "    The MPI communicator.                                           \n"
        )
        .def("get_voltage", (&pycap::get_voltage),
            "Measure the voltage across the device.                         \n"
            "                                                               \n"
            "Returns                                                        \n"
            "-------                                                        \n"
            "float                                                          \n"
            "    The voltage across the device in units of volts.           \n"
            ,
            boost::python::args("self") )
        .def("get_current", (&pycap::get_current),
            "Measure the electrical current that flows through the device.  \n"
            "                                                               \n"
            "Returns                                                        \n"
            "-------                                                        \n"
            "float                                                          \n"
            "    The current in units of amperes.                           \n"
            ,
            boost::python::args("self") )
        .def("evolve_one_time_step_constant_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_current),
            "Impose the electrical current and evolve in time.              \n"
            "                                                               \n"
            "Parameters                                                     \n"
            "----------                                                     \n"
            "time_step : float                                              \n"
            "    The time step in units of seconds.                         \n"
            "current   : float                                              \n"
            "    The electrical current in units of amperes.                \n"
            ,
            boost::python::args("self", "time_step", "current") )
        .def("evolve_one_time_step_constant_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_voltage),
            "Impose the voltage across the device and evolve in time.       \n"
            "                                                               \n"
            "Parameters                                                     \n"
            "----------                                                     \n"
            "time_step : float                                              \n"
            "    The time step in units of seconds.                         \n"
            "voltage   : float                                              \n"
            "    The voltage across the device in units of volts.           \n"
            ,
            boost::python::args("self", "time_step", "voltage") )
        .def("evolve_one_time_step_constant_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_power  ),
            "Impose the power and evolve in time.                           \n"
            "                                                               \n"
            "Parameters                                                     \n"
            "----------                                                     \n"
            "time_step : float                                              \n"
            "    The time step in units of seconds.                         \n"
            "power     : float                                              \n"
            "    The power in units of watts.                               \n"
            ,
            boost::python::args("self", "time_step", "power"  ) )
        .def("evolve_one_time_step_constant_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_load   ),
            "Impose the load and evolve in time.                            \n"
            "                                                               \n"
            "Parameters                                                     \n"
            "----------                                                     \n"
            "time_step : float                                              \n"
            "    The time step in units of seconds.                         \n"
            "load      : float                                              \n"
            "    The load in units of ohms.                                 \n"
            ,
            boost::python::args("self", "time_step", "power"  ) )
        .def("evolve_one_time_step_linear_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_current), boost::python::args("self", "time_step", "current") )
        .def("evolve_one_time_step_linear_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_voltage), boost::python::args("self", "time_step", "voltage") )
        .def("evolve_one_time_step_linear_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_power  ), boost::python::args("self", "time_step", "power"  ) )
        .def("evolve_one_time_step_linear_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_load   ), boost::python::args("self", "time_step", "load"   ) )
        .def("compute_equivalent_circuit", &pycap::compute_equivalent_circuit,
            "Compute the equivalent circuit to a supercapacitor.            \n"
            "                                                               \n"
            "Parameters                                                     \n"
            "----------                                                     \n"
            "ptree : pycap.PropertyTree                                     \n"
            "    The tree to build a supercapacitor.                        \n"
            "                                                               \n"
            "Returns                                                        \n"
            "-------                                                        \n"
            "pycap.PropertyTree                                             \n"
            "    The tree to build the equivalent circuit.                  \n"
            "                                                               \n"
            "Examples                                                       \n"
            "--------                                                       \n"
            ">>> from pycap import PropertyTree, EnergyStorageDevice        \n"
            ">>> super_capacitor_ptree=PropertyTree()                       \n"
            ">>> super_capacitor_ptree.parse_info('super_capacitor.info')   \n"
            ">>> super_capacitor = EnergyStorageDevice(super_capacitor_ptree)\n"
            ">>> equivalent_circuit_ptree = EnergyStorageDevice.compute_equivalent_circuit(super_capacitor_ptree)\n"
            ">>> equivalent_circuit = EnergyStorageDevice(equivalent_circuit_ptree)\n"
            ,
            boost::python::args("ptree") )
        .staticmethod("compute_equivalent_circuit")
//        .def_pickle(pycap::serializable_class_pickle_support<cap::EnergyStorageDevice>())
        ;

    boost::python::class_<
        boost::property_tree::ptree,
        std::shared_ptr<boost::property_tree::ptree> >(
        "PropertyTree",
        "Wrappers for Boost.PropertyTree\n"
        "\n"
        "Examples\n"
        "--------\n"
        ">>> ptree = PropertyTree()\n"
        ">>> ptree.put_double('pi', 3.14)\n"
        ">>> ptree.get_double('pi')\n"
        "3.14\n"
        ">>> ptree.get_double_with_default('sqrt2', 1,41)\n"
        "1.41\n"
        "\n"
        "Raises\n"
        "------\n"
        "RuntimeError: No such node (<path>)\n"
        "    Error indicating that specified <path> does not exist.\n"
        "RuntimeError: conversion of data to type \"<type>\" failed\n"
        "    Error indicating that translation from or to <type> has\n"
        "    failed.\n"
        )
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
        .def("parse_info"                   , &pycap::parse_info                   , "Read the input file at INFO format and populate the PropertyTree.", boost::python::args("self", "filename") )
        .def("get_child"                    , &pycap::get_child                    , "Get the child at the given path, or throw ptree_bad_path."        , boost::python::args("self", "path") )
        .def_pickle(pycap::serializable_class_pickle_support<boost::property_tree::ptree>())
        ;
}
