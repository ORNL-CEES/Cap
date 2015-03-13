#ifndef CAP_ENERGY_STORAGE_DEVICE_H
#define CAP_ENERGY_STORAGE_DEVICE_H

#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace cap {

class Parameters
{
public:
    Parameters(std::shared_ptr<boost::property_tree::ptree const> d) : database(d) { }

    std::shared_ptr<boost::property_tree::ptree const> database;
};

class EnergyStorageDevice
{
public:
    EnergyStorageDevice() { }
    virtual ~EnergyStorageDevice() { }
    virtual void   print_data(std::ostream & os) const = 0;
    virtual double get_voltage() const = 0;
    virtual double get_current() const = 0;
    virtual void   reset_voltage(double const voltage) = 0;
    virtual void   reset_current(double const current) = 0;
    virtual void   evolve_one_time_step_constant_current(double const time_step, double const constant_current) = 0;
    virtual void   evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage) = 0;
    virtual void   evolve_one_time_step_constant_power  (double const time_step, double const constant_power  ) = 0;
};



std::shared_ptr<cap::EnergyStorageDevice>
buildEnergyStorageDevice(std::shared_ptr<cap::Parameters const> params);

}

#endif // CAP_ENERGY_STORAGE_DEVICE_H
