/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_ENERGY_STORAGE_DEVICE_H
#define CAP_ENERGY_STORAGE_DEVICE_H

#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/access.hpp>
#include <boost/mpi/communicator.hpp>
#include <memory>
#include <map>

namespace cap
{

class EnergyStorageDeviceBuilder;
class EnergyStorageDeviceInspector;

/**
 * This class is an abstract representation of an energy storage device. It can
 * evolve in time at various operating conditions and return the voltage drop
 * across itself and the electrical current that flows through it.
 */
class EnergyStorageDevice
{
public:
  EnergyStorageDevice(boost::mpi::communicator const &communicator);

  virtual ~EnergyStorageDevice();

  virtual void get_voltage(double &voltage) const = 0;

  virtual void get_current(double &current) const = 0;

  /**
   * Advance the time by @p time_step seconds. The current is constant during
   * time step and its value is @p current amperes.
   */
  virtual void evolve_one_time_step_constant_current(double const time_step,
                                                     double const current) = 0;

  /**
   * Advance the time by @p time_step seconds. The voltage is constant during
   * time step and its value is @p voltage volts.
   */
  virtual void evolve_one_time_step_constant_voltage(double const time_step,
                                                     double const voltage) = 0;

  /**
   * Advance the time by @p time_step seconds. The power is constant during
   * time step and its value is @p power watts.
   */
  virtual void evolve_one_time_step_constant_power(double const time_step,
                                                   double const power) = 0;

  /**
   * Advance the time by @p time_step seconds. The load is constant during
   * time step and its value is @p load ohms.
   */
  virtual void evolve_one_time_step_constant_load(double const time_step,
                                                  double const load) = 0;

  /**
   * Advance the time by @p time_step seconds. The current is changing linearly
   * during time step and its value at the end of the time step is
   * @p current amperes.
   */
  virtual void evolve_one_time_step_linear_current(double const time_step,
                                                   double const current) = 0;

  /**
   * Advance the time by @p time_step seconds. The voltage is changing linearly
   * during time step and its value at the end of the time step is
   * @p voltage volts.
   */
  virtual void evolve_one_time_step_linear_voltage(double const time_step,
                                                   double const voltage) = 0;

  /**
   * Advance the time by @p time_step seconds. The power is changing linearly
   * during time step and its value at the end of the time step is
   * @p power watts.
   */
  virtual void evolve_one_time_step_linear_power(double const time_step,
                                                 double const power) = 0;

  /**
   * Advance the time by @p time_step seconds. The load is changing linearly
   * during time step and its value at the end of the time step is
   * @p load ohms.
   */
  virtual void evolve_one_time_step_linear_load(double const time_step,
                                                double const load) = 0;

  /**
   * Factory function that creates an EnergyStorageDevice object.
   */
  static std::unique_ptr<EnergyStorageDevice>
  build(boost::mpi::communicator const &comm,
        boost::property_tree::ptree const &ptree);

  /**
   * Visitor design pattern of the EnergyStorageDevice object. This allows to
   * define a new operation without changing the classes of the elements on
   * which it operates.
   */
  virtual void inspect(EnergyStorageDeviceInspector *inspector) = 0;

protected:
  boost::mpi::communicator _communicator;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    // nothing to do
    std::ignore = ar;
    std::ignore = version;
  }
  friend EnergyStorageDeviceBuilder;
  static std::map<std::string, EnergyStorageDeviceBuilder *> _builders;
};

/**
 * Helper class for the factory function.
 */
class EnergyStorageDeviceBuilder
{
public:
  virtual ~EnergyStorageDeviceBuilder();
  virtual std::unique_ptr<EnergyStorageDevice>
  build(boost::mpi::communicator const &comm,
        boost::property_tree::ptree const &ptree) = 0;

  /**
   * This function registers the energy device storages so that they can be
   * called by the factory function.
   */
  static void
  register_energy_storage_device(std::string const &type,
                                 EnergyStorageDeviceBuilder *builder);
};

/**
 * Abstract class used for the visitor design pattern.
 */
class EnergyStorageDeviceInspector
{
public:
  virtual ~EnergyStorageDeviceInspector();
  virtual void inspect(EnergyStorageDevice *device) = 0;
};

#define REGISTER_ENERGY_STORAGE_DEVICE(T)                                      \
  class T##Builder : public EnergyStorageDeviceBuilder                         \
  {                                                                            \
  public:                                                                      \
    T##Builder() { register_energy_storage_device(#T, this); }                 \
    virtual std::unique_ptr<EnergyStorageDevice>                               \
    build(boost::mpi::communicator const &comm,                                \
          boost::property_tree::ptree const &ptree) override                   \
    {                                                                          \
      return std::unique_ptr<T>(new T(comm, ptree));                           \
    }                                                                          \
  };                                                                           \
  static T##Builder global_##T##Builder;
} // end namespace cap

#endif // CAP_ENERGY_STORAGE_DEVICE_H
