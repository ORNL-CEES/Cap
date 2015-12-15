#ifndef CAP_RESISTOR_CAPACITOR_H
#define CAP_RESISTOR_CAPACITOR_H

#include <cap/energy_storage_device.h>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <string>

namespace cap {

class SeriesRC : public EnergyStorageDevice
{
public:
    SeriesRC(boost::mpi::communicator const & comm, boost::property_tree::ptree const & ptree);
    void inspect(EnergyStorageDeviceInspector * inspector) override;
    void print_data(std::ostream & os) const override;
    void reset_voltage(double const voltage) override;
    void reset_current(double const current) override;
    void evolve_one_time_step_constant_current(double const delta_t, double const constant_current) override;
    void evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage) override;
    void evolve_one_time_step_constant_power  (double const delta_t, double const constant_power  ) override;
    void evolve_one_time_step_constant_load   (double const delta_t, double const constant_load   ) override;
    void evolve_one_time_step_changing_current(double const delta_t, double const changing_current) override;
    void evolve_one_time_step_changing_voltage(double const delta_t, double const changing_voltage) override;
    inline void get_voltage(double & voltage) const override { voltage = U; }
    inline void get_current(double & current) const override { current = I; }
    // TODO:
    std::size_t evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method = "NEWTON");
    void reset(double const capacitor_voltage);
    
    double R;
    double C;
    double U_C;
    double U;
    double I;
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<cap::EnergyStorageDevice>(*this);
        ar & R & C & U_C & U & I;
        std::ignore = version;
    }
};

class ParallelRC : public EnergyStorageDevice
{
public:
    ParallelRC(boost::mpi::communicator const & comm, boost::property_tree::ptree const & ptree);
    void inspect(EnergyStorageDeviceInspector * inspector) override;
    void print_data(std::ostream & os) const override;
    void reset_voltage(double const voltage) override;
    void reset_current(double const current) override;
    void evolve_one_time_step_constant_current(double const delta_t, double const constant_current) override;
    void evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage) override;
    void evolve_one_time_step_constant_power  (double const delta_t, double const constant_power  ) override;
    void evolve_one_time_step_constant_load   (double const delta_t, double const constant_load   ) override;
    void evolve_one_time_step_changing_current(double const delta_t, double const changing_current) override;
    void evolve_one_time_step_changing_voltage(double const delta_t, double const constant_voltage) override;
    inline void get_voltage(double & voltage) const override { voltage = U; }
    inline void get_current(double & current) const override { current = I; }
    // TODO:
    std::size_t evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method = "NEWTON");
    void reset(double const capacitor_voltage);

    double R_parallel;
    double C;
    double U_C;
    double U;
    double I;
    double R_series;
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<cap::EnergyStorageDevice>(*this);
        ar & R_parallel & C & U_C & U & I & R_series;
        std::ignore = version;
    }
};

} // end namespace cap

namespace boost { namespace serialization {
template<class Archive>
inline void save_construct_data(
    Archive & ar, const cap::SeriesRC * rc, const unsigned int file_version
){
    std::ignore = file_version;
    ar << rc->R << rc->C << rc->U_C << rc->U << rc->I;
}

template<class Archive>
inline void load_construct_data(
    Archive & ar, cap::SeriesRC * rc, const unsigned int file_version
){
    std::ignore = file_version;
    double R, C, U_C, U, I;
    ar >> R >> C >> U_C >> U >> I;
    boost::property_tree::ptree ptree;
    ptree.put("series_resistance", R);
    ptree.put("capacitance"      , C);
    ::new(rc)cap::SeriesRC(boost::mpi::communicator(), ptree);
}

template<class Archive>
inline void save_construct_data(
    Archive & ar, const cap::ParallelRC * rc, const unsigned int file_version
){
    std::ignore = file_version;
    ar << rc->R_parallel << rc->C << rc->U_C << rc->U << rc->I << rc->R_series;
}

template<class Archive>
inline void load_construct_data(
    Archive & ar, cap::ParallelRC * rc, const unsigned int file_version
){
    std::ignore = file_version;
    double R_parallel, C, U_C, U, I, R_series;
    ar >> R_parallel >> C >> U_C >> U >> I >> R_series;
    boost::property_tree::ptree ptree;
    ptree.put("series_resistance", R_series);
    ptree.put("parallel_resistance", R_parallel);
    ptree.put("capacitance", C);
    ::new(rc)cap::ParallelRC(boost::mpi::communicator(), ptree);
}
}} // namespace ...

#endif // CAP_RESISTOR_CAPACITOR_H
