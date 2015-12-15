#ifndef CAP_ENERGY_STORAGE_DEVICE_H
#define CAP_ENERGY_STORAGE_DEVICE_H

#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/access.hpp>
#include <boost/mpi/communicator.hpp>
#include <memory>

namespace cap {

class EnergyStorageDevice
{
public:
    EnergyStorageDevice(boost::mpi::communicator const & communicator);
    virtual ~EnergyStorageDevice();
    virtual void print_data(std::ostream & os) const = 0;
    virtual void get_voltage(double & voltage) const = 0;
    virtual void get_current(double & current) const = 0;
    virtual void reset_voltage(double const voltage) = 0;
    virtual void reset_current(double const current) = 0;
    virtual void evolve_one_time_step_constant_current(double const time_step, double const constant_current) = 0;
    virtual void evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage) = 0;
    virtual void evolve_one_time_step_constant_power  (double const time_step, double const constant_power  ) = 0;
    virtual void evolve_one_time_step_constant_load   (double const time_step, double const constant_load   ) = 0;
    virtual void evolve_one_time_step_changing_current(double const time_step, double const changing_current);
    virtual void evolve_one_time_step_changing_voltage(double const time_step, double const changing_voltage);
    virtual void evolve_one_time_step_changing_power  (double const time_step, double const changing_power  );
    virtual void evolve_one_time_step_changing_load   (double const time_step, double const changing_load   );
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // nothing to do
        std::ignore = ar;
        std::ignore = version;
    }
    boost::mpi::communicator communicator_;
};

std::shared_ptr<EnergyStorageDevice>
buildEnergyStorageDevice(boost::mpi::communicator const & communicator,
                         boost::property_tree::ptree const & ptree);

}

#endif // CAP_ENERGY_STORAGE_DEVICE_H
