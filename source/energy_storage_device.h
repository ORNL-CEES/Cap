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

class EnergyStorageDevice
{
public:
  EnergyStorageDevice(boost::mpi::communicator const &communicator);
  virtual ~EnergyStorageDevice();
  // DEPRECATED ///////////////////////////////////////
  virtual void print_data(std::ostream &os) const = 0;
  virtual void reset_voltage(double const voltage) = 0;
  virtual void reset_current(double const current) = 0;
  /////////////////////////////////////////////////////
  virtual void get_voltage(double &voltage) const = 0;
  virtual void get_current(double &current) const = 0;
  virtual void
  evolve_one_time_step_constant_current(double const time_step,
                                        double const constant_current) = 0;
  virtual void
  evolve_one_time_step_constant_voltage(double const time_step,
                                        double const constant_voltage) = 0;
  virtual void
  evolve_one_time_step_constant_power(double const time_step,
                                      double const constant_power) = 0;
  virtual void
  evolve_one_time_step_constant_load(double const time_step,
                                     double const constant_load) = 0;
  virtual void
  evolve_one_time_step_changing_current(double const time_step,
                                        double const changing_current);
  virtual void
  evolve_one_time_step_changing_voltage(double const time_step,
                                        double const changing_voltage);
  virtual void evolve_one_time_step_changing_power(double const time_step,
                                                   double const changing_power);
  virtual void evolve_one_time_step_changing_load(double const time_step,
                                                  double const changing_load);
  static std::unique_ptr<EnergyStorageDevice>
  build(boost::mpi::communicator const &comm,
        boost::property_tree::ptree const &ptree);
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

class EnergyStorageDeviceBuilder
{
public:
  virtual ~EnergyStorageDeviceBuilder();
  virtual std::unique_ptr<EnergyStorageDevice>
  build(boost::mpi::communicator const &comm,
        boost::property_tree::ptree const &ptree) = 0;
  static void
  register_energy_storage_device(std::string const &type,
                                 EnergyStorageDeviceBuilder *builder);
};

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

// DEPRECATED /////////////////////////////////////////////////////////
std::shared_ptr<EnergyStorageDevice>
buildEnergyStorageDevice(boost::mpi::communicator const &communicator,
                         boost::property_tree::ptree const &ptree);
///////////////////////////////////////////////////////////////////////

} // end namespace cap

#endif // CAP_ENERGY_STORAGE_DEVICE_H
