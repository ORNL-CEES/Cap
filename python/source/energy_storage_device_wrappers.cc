#include <pycap/energy_storage_device_wrappers.h>
#include <cap/equivalent_circuit.h>
#include <mpi4py/mpi4py.h>

namespace pycap {

void EnergyStorageDeviceWrap::evolve_one_time_step_constant_current(double const time_step, double const current)
{
    this->get_override("evolve_one_time_step_constant_current")(time_step, current);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_constant_voltage(double const time_step, double const voltage)
{
    this->get_override("evolve_one_time_step_constant_voltage")(time_step, voltage);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_constant_power  (double const time_step, double const power  )
{
    this->get_override("evolve_one_time_step_constant_power"  )(time_step, power  );
}
void EnergyStorageDeviceWrap::evolve_one_time_step_constant_load   (double const time_step, double const load   )
{
    this->get_override("evolve_one_time_step_constant_load"   )(time_step, load   );
}
void EnergyStorageDeviceWrap::evolve_one_time_step_linear_current  (double const time_step, double const current)
{
    this->get_override("evolve_one_time_step_new_current"     )(time_step, current);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_linear_voltage  (double const time_step, double const voltage)
{
    this->get_override("evolve_one_time_step_linear_voltage"  )(time_step, voltage);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_linear_power    (double const time_step, double const power  )
{
    this->get_override("evolve_one_time_step_linear_power"    )(time_step, power  );
}
void EnergyStorageDeviceWrap::evolve_one_time_step_linear_load     (double const time_step, double const load   )
{
    this->get_override("evolve_one_time_step_linear_load"     )(time_step, load   );
}

double get_current(cap::EnergyStorageDevice const & dev)
{
    double current;
    dev.get_current(current);
    return current;
}

double get_voltage(cap::EnergyStorageDevice const & dev)
{
    double voltage;
    dev.get_voltage(voltage);
    return voltage;
}

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & py_comm, boost::python::object & py_ptree)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(py_ptree);
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
    if (comm_p == NULL) boost::python::throw_error_already_set();
    boost::mpi::communicator comm(*comm_p, boost::mpi::comm_attach);
    return cap::EnergyStorageDevice::build(comm, ptree);
}

std::shared_ptr<boost::property_tree::ptree>
compute_equivalent_circuit(boost::python::object & python_object)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(python_object);
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(ptree);
    std::shared_ptr<boost::property_tree::ptree> equivalent_circuit_database =
        std::make_shared<boost::property_tree::ptree>();
    cap::compute_equivalent_circuit(device_database, equivalent_circuit_database);
    return equivalent_circuit_database;
}
} // end namespace pycap

