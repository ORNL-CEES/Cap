#ifndef ENERGY_STORAGE_DEVICE_WRAPPERS_H
#define ENERGY_STORAGE_DEVICE_WRAPPERS_H

#include <cap/energy_storage_device.h>
#include <boost/python/object.hpp>
#include <boost/python/wrapper.hpp>

namespace pycap {

struct EnergyStorageDeviceWrap : cap::EnergyStorageDevice, boost::python::wrapper<cap::EnergyStorageDevice>
{
    void evolve_one_time_step_constant_current(double const time_step, double const current);
    void evolve_one_time_step_constant_voltage(double const time_step, double const voltage);
    void evolve_one_time_step_constant_power  (double const time_step, double const power  );
    void evolve_one_time_step_constant_load   (double const time_step, double const load   );
    void evolve_one_time_step_linear_current  (double const time_step, double const current);
    void evolve_one_time_step_linear_voltage  (double const time_step, double const voltage);
    void evolve_one_time_step_linear_power    (double const time_step, double const power  );
    void evolve_one_time_step_linear_load     (double const time_step, double const load   );
};

double get_current(cap::EnergyStorageDevice const & device);
double get_voltage(cap::EnergyStorageDevice const & device);

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & py_ptree,
                            boost::python::object & py_comm);

boost::property_tree::ptree compute_equivalent_circuit(boost::python::object & python_object);

} // end namespace pycap

#endif
