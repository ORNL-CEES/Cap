/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_RESISTOR_CAPACITOR_H
#define CAP_RESISTOR_CAPACITOR_H

#include <cap/energy_storage_device.h>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <string>

namespace cap
{

class SeriesRC : public EnergyStorageDevice
{
public:
  SeriesRC(boost::property_tree::ptree const &ptree,
           boost::mpi::communicator const &comm);

  void inspect(EnergyStorageDeviceInspector *inspector) override;

  void evolve_one_time_step_constant_current(double const delta_t,
                                             double const current) override;

  void evolve_one_time_step_constant_voltage(double const delta_t,
                                             double const voltage) override;

  void evolve_one_time_step_constant_power(double const delta_t,
                                           double const power) override;

  void evolve_one_time_step_constant_load(double const delta_t,
                                          double const load) override;

  void evolve_one_time_step_linear_current(double const delta_t,
                                           double const current) override;

  void evolve_one_time_step_linear_voltage(double const delta_t,
                                           double const voltage) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_power(double const delta_t,
                                         double const power) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_load(double const delta_t,
                                        double const load) override;

  void get_voltage(double &voltage) const override { voltage = U; }

  void get_current(double &current) const override { current = I; }

  /**
   * This function advance the time by @p delta_t seconds. The power is
   * constant during the time step and its value is @p power. This
   * function internally solves a non-linear problem and the third parameter @p
   * method allows to choose between FIXED_POINT (Picard iteration) and NEWTON.
   * This function returns the number of non-linear iterations performed to
   * reach convergence.
   */
  std::size_t
  evolve_one_time_step_constant_power(double const delta_t, double const power,
                                      std::string const &method = "NEWTON");

  /**
   * Save the current state of energy device in a file.
   */
  void save(const std::string &filename) const override;

  /**
   * Load an energy device from a state saved in a file.
   */
  void load(const std::string &filename) override;

  // TODO: make these variables private
  double R;
  double C;
  double U_C;
  double U;
  double I;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &boost::serialization::base_object<cap::EnergyStorageDevice>(*this);
    ar &R &C &U_C &U &I;
    std::ignore = version;
  }

  boost::mpi::communicator _comm;
};

class ParallelRC : public EnergyStorageDevice
{
public:
  ParallelRC(boost::property_tree::ptree const &ptree,
             boost::mpi::communicator const &comm);

  void inspect(EnergyStorageDeviceInspector *inspector) override;

  void evolve_one_time_step_constant_current(double const delta_t,
                                             double const current) override;

  void evolve_one_time_step_constant_voltage(double const delta_t,
                                             double const voltage) override;

  void evolve_one_time_step_constant_power(double const delta_t,
                                           double const power) override;

  void evolve_one_time_step_constant_load(double const delta_t,
                                          double const load) override;

  void evolve_one_time_step_linear_current(double const delta_t,
                                           double const current) override;

  void evolve_one_time_step_linear_voltage(double const delta_t,
                                           double const voltage) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_power(double const delta_t,
                                         double const power) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_load(double const delta_t,
                                        double const load) override;

  inline void get_voltage(double &voltage) const override { voltage = U; }

  inline void get_current(double &current) const override { current = I; }

  /**
   * This function advance the time by @p delta_t seconds. The power is
   * constant during the time step and its value is @p power. This
   * function internally solves a non-linear problem and the third parameter @p
   * method allows to choose between FIXED_POINT (Picard iteration) and NEWTON.
   * This function returns the number of non-linear iterations performed to
   * reach convergence.
   */
  std::size_t
  evolve_one_time_step_constant_power(double const delta_t, double const power,
                                      std::string const &method = "NEWTON");

  /**
   * Save the current state of energy device in a file.
   */
  void save(const std::string &filename) const override;

  /**
   * Load an energy device from a state saved in a file.
   */
  void load(const std::string &filename) override;

  // TODO: make these variables private
  double R_series;
  double R_parallel;
  double C;
  double U_C;
  double U;
  double I;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar &boost::serialization::base_object<cap::EnergyStorageDevice>(*this);
    ar &R_parallel &C &U_C &U &I &R_series;
    std::ignore = version;
  }

  boost::mpi::communicator _comm;
};

} // end namespace cap

namespace boost
{
namespace serialization
{
/**
 * Save the current state of a cap::SeriesRC object.
 */
template <class Archive>
inline void save_construct_data(Archive &ar, const cap::SeriesRC *rc,
                                const unsigned int file_version)
{
  std::ignore = file_version;
  ar << rc->R << rc->C << rc->U_C << rc->U << rc->I;
}

/**
 * Build a cap::SeriesRC object from saved data.
 */
template <class Archive>
inline void load_construct_data(Archive &ar, cap::SeriesRC *rc,
                                const unsigned int file_version)
{
  std::ignore = file_version;
  double R, C, U_C, U, I;
  ar >> R >> C >> U_C >> U >> I;
  boost::property_tree::ptree ptree;
  ptree.put("series_resistance", R);
  ptree.put("capacitance", C);
  ::new (rc) cap::SeriesRC(ptree, boost::mpi::communicator());
}

/**
 * Save the current state of a cap::ParallelRC object.
 */
template <class Archive>
inline void save_construct_data(Archive &ar, const cap::ParallelRC *rc,
                                const unsigned int file_version)
{
  std::ignore = file_version;
  ar << rc->R_parallel << rc->C << rc->U_C << rc->U << rc->I << rc->R_series;
}

/**
 * Build a cap::ParallelRC object from saved data.
 */
template <class Archive>
inline void load_construct_data(Archive &ar, cap::ParallelRC *rc,
                                const unsigned int file_version)
{
  std::ignore = file_version;
  double R_parallel, C, U_C, U, I, R_series;
  ar >> R_parallel >> C >> U_C >> U >> I >> R_series;
  boost::property_tree::ptree ptree;
  ptree.put("series_resistance", R_series);
  ptree.put("parallel_resistance", R_parallel);
  ptree.put("capacitance", C);
  ::new (rc) cap::ParallelRC(ptree, boost::mpi::communicator());
}
}
} // namespace ...

#endif // CAP_RESISTOR_CAPACITOR_H
