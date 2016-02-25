#include <pycap/energy_storage_device_wrappers.h>
#include <boost/python.hpp>
#include <mpi4py/mpi4py.h>

namespace pycap
{

char const energy_storage_device_docstring[] =
  "Wrappers for Cap.EnergyStorageDevice                                     \n"
  "                                                                         \n"
  "Examples                                                                 \n"
  "--------                                                                 \n"
  "                                                                         \n"
  ">>> ptree = PropertyTree()                                               \n"
  ">>> ptree.parse_info('device.info')                                      \n"
  ">>> device = EnergyStorageDevice(ptree)                                  \n"
  ">>> delta_t, U = 1.0, 2.1 # units of seconds and volts                   \n"
  ">>> device.evolve_one_time_step_constant_voltage(delta_t,U)              \n"
  ">>> I = device.get_current() # amperes                                   \n"
  "                                                                         \n"
  ;

char const get_voltage_docstring[] =
  "Measure the voltage across the device.                                   \n"
  "                                                                         \n"
  "Returns                                                                  \n"
  "-------                                                                  \n"
  "float                                                                    \n"
  "    The voltage across the device in volts.                              \n"
  ;

char const get_current_docstring[] =
  "Measure the electrical current that flows through the device.            \n"
  "                                                                         \n"
  "Returns                                                                  \n"
  "-------                                                                  \n"
  "float                                                                    \n"
  "    The electrical current in amperes.                                   \n"
  ;

char const evolve_one_time_step_constant_current_docstring[] =
  "Impose the electrical current and evolve in time.                        \n"
  "                                                                         \n"
  "Parameters                                                               \n"
  "----------                                                               \n"
  "time_step : float                                                        \n"
  "    The time step in seconds.                                            \n"
  "current : float                                                          \n"
  "    The electrical current in amperes.                                   \n"
  ;

char const evolve_one_time_step_constant_voltage_docstring[] =
  "Impose the voltage across the device and evolve in time.                 \n"
  "                                                                         \n"
  "Parameters                                                               \n"
  "----------                                                               \n"
  "time_step : float                                                        \n"
  "    The time step in seconds.                                            \n"
  "voltage : float                                                          \n"
  "    The voltage across the device in volts.                              \n"
  ;

char const evolve_one_time_step_constant_power_docstring[] =
  "Impose the power and evolve in time.                                     \n"
  "                                                                         \n"
  "Parameters                                                               \n"
  "----------                                                               \n"
  "time_step : float                                                        \n"
  "    The time step in seconds.                                            \n"
  "power : float                                                            \n"
  "    The power in watts.                                                  \n"
  ;

char const evolve_one_time_step_constant_load_docstring[] =
  "Impose the load and evolve in time.                                      \n"
  "                                                                         \n"
  "Parameters                                                               \n"
  "----------                                                               \n"
  "time_step : float                                                        \n"
  "    The time step in seconds.                                            \n"
  "load : float                                                             \n"
  "    The load in ohms.                                                    \n"
  ;

char const compute_equivalent_circuit_docstring[] =
  "Compute the equivalent circuit to a supercapacitor.                      \n"
  "                                                                         \n"
  "Parameters                                                               \n"
  "----------                                                               \n"
  "ptree : pycap.PropertyTree                                               \n"
  "    The tree to build a supercapacitor.                                  \n"
  "                                                                         \n"
  "Returns                                                                  \n"
  "-------                                                                  \n"
  "pycap.PropertyTree                                                       \n"
  "    The tree to build the equivalent circuit.                            \n"
  "                                                                         \n"
  "Examples                                                                 \n"
  "--------                                                                 \n"
  ">>> from pycap import PropertyTree, EnergyStorageDevice                  \n"
  ">>> super_capacitor_ptree = PropertyTree()                               \n"
  ">>> super_capacitor_ptree.parse_info('super_capacitor.info')             \n"
  ">>> super_capacitor = EnergyStorageDevice(super_capacitor_ptree)         \n"
  ">>> equivalent_circuit_ptree = EnergyStorageDevice.compute_equivalent_circuit(super_capacitor_ptree)\n"
  ">>> equivalent_circuit = EnergyStorageDevice(equivalent_circuit_ptree)   \n"
  ;

void export_energy_storage_device()
{
  boost::python::class_<EnergyStorageDeviceWrap,
                        std::shared_ptr<EnergyStorageDeviceWrap>,
                        boost::noncopyable> (
    "EnergyStorageDevice",
    energy_storage_device_docstring,
    boost::python::no_init)
    .def("build",
    boost::python::make_constructor(&build_energy_storage_device,
    boost::python::default_call_policies(),
    boost::python::args("ptree", "comm")),
    "                                                                       \n"
    "Parameters                                                             \n"
    "----------                                                             \n"
    "ptree : pycap.PropertyTree                                             \n"
    "    The appropriate property tree to create a device                   \n"
    "    from the factory.                                                  \n"
    "comm  : mpi4py.MPI.Comm                                                \n"
    "    The MPI communicator.                                              \n"
    )
    .staticmethod("build")
    .def("get_voltage", &get_voltage, get_voltage_docstring,
         boost::python::args("self") )
    .def("get_current", &get_current, get_current_docstring,
         boost::python::args("self") )
    .def("evolve_one_time_step_constant_current",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_current),
         evolve_one_time_step_constant_current_docstring,
         boost::python::args("self", "time_step", "current") )
    .def("evolve_one_time_step_constant_voltage",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_voltage),
         evolve_one_time_step_constant_voltage_docstring,
         boost::python::args("self", "time_step", "voltage") )
    .def("evolve_one_time_step_constant_power",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_power),
         evolve_one_time_step_constant_power_docstring,
         boost::python::args("self", "time_step", "power") )
    .def("evolve_one_time_step_constant_load",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_load),
         evolve_one_time_step_constant_load_docstring,
         boost::python::args("self", "time_step", "load") )
    .def("evolve_one_time_step_linear_current",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_current),
         boost::python::args("self", "time_step", "current") )
    .def("evolve_one_time_step_linear_voltage",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_voltage),
         boost::python::args("self", "time_step", "voltage") )
    .def("evolve_one_time_step_linear_power",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_power),
         boost::python::args("self", "time_step", "power") )
    .def("evolve_one_time_step_linear_load",
         boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_linear_load),
         boost::python::args("self", "time_step", "load") )
    .def("compute_equivalent_circuit", &compute_equivalent_circuit,
         compute_equivalent_circuit_docstring,
         boost::python::args("ptree") )
    .staticmethod("compute_equivalent_circuit")
//        .def_pickle(pycap::serializable_class_pickle_support<cap::EnergyStorageDevice>())
        ;
}

} // end namespace pycap
