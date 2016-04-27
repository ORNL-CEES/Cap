#ifndef ENERGY_STORAGE_DEVICE_WRAPPERS_H
#define ENERGY_STORAGE_DEVICE_WRAPPERS_H

#include <cap/energy_storage_device.h>
#include <boost/python/object.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/dict.hpp>

namespace pycap {

double get_current(cap::EnergyStorageDevice const & device);
double get_voltage(cap::EnergyStorageDevice const & device);
// TODO: may want const reference here
boost::python::dict inspect(cap::EnergyStorageDevice & device);

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & py_ptree,
                            boost::python::object & py_comm);

} // end namespace pycap

#endif
