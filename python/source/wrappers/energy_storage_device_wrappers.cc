/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <pycap/energy_storage_device_wrappers.h>
#include <cap/default_inspector.h>
#include <mpi4py/mpi4py.h>

namespace pycap {

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

boost::python::dict inspect(cap::EnergyStorageDevice & dev)
{
    cap::DefaultInspector inspector;
    dev.inspect(&inspector);
    boost::python::dict data;
    for (auto x : inspector.get_data())
        data[x.first] = x.second;
    return data;
}

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & py_ptree,
                            boost::python::object & py_comm)
{
    if (import_mpi4py() < 0) throw std::runtime_error("Failed to import mpi4py");
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(py_ptree);
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
    if (comm_p == nullptr) boost::python::throw_error_already_set();
    boost::mpi::communicator comm(*comm_p, boost::mpi::comm_attach);
    return cap::EnergyStorageDevice::build(ptree, comm);
}

} // end namespace pycap

