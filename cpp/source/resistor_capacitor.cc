#include <cap/resistor_capacitor.h>
#include <boost/format.hpp>
#include <cmath>
#include <stdexcept>

namespace cap
{

REGISTER_ENERGY_STORAGE_DEVICE(SeriesRC)
REGISTER_ENERGY_STORAGE_DEVICE(ParallelRC)

void ParallelRC::inspect(EnergyStorageDeviceInspector *inspector)
{
  inspector->inspect(this);
}

void SeriesRC::inspect(EnergyStorageDeviceInspector *inspector)
{
  inspector->inspect(this);
}

SeriesRC::SeriesRC(boost::mpi::communicator const &comm,
                   boost::property_tree::ptree const &ptree)
    : EnergyStorageDevice(comm)
{
  R = ptree.get<double>("series_resistance");
  C = ptree.get<double>("capacitance");
  U = ptree.get<double>("initial_voltage", 0.0);
  this->reset_voltage(U);
}

void SeriesRC::print_data(std::ostream &os) const
{
  os << boost::format("  %10.7f  %10.7f  %10.7f  \n") % I % U % U_C;
}

void SeriesRC::reset_voltage(double const voltage)
{
  U_C = voltage;
  U   = U_C;
  I   = 0.0;
}

void SeriesRC::reset_current(double const current)
{
  I = current;
  U = U_C + R * I;
}

void SeriesRC::reset(double const capacitor_voltage)
{
  U_C = capacitor_voltage;
  U   = U_C + R * I;
}

void SeriesRC::evolve_one_time_step_constant_load(double const delta_t,
                                                  double const load)
{
  U_C *= std::exp(-delta_t / ((R + load) * C));
  I = -U_C / (R + load);
  U = U_C + R * I;
}

void SeriesRC::evolve_one_time_step_constant_power(double const delta_t,
                                                   double const power)
{
  evolve_one_time_step_constant_power(delta_t, power, "NEWTON");
}

void SeriesRC::evolve_one_time_step_constant_current(double const delta_t,
                                                     double const current)
{
  U_C += current * delta_t / C;
  I = current;
  U = R * I + U_C;
}

void SeriesRC::evolve_one_time_step_linear_current(double const delta_t,
                                                   double const current)
{
  U_C += I * delta_t / C;
  U_C += (current - I) * 0.5 * delta_t / C;
  I = current;
  U = R * I + U_C;
}

void SeriesRC::evolve_one_time_step_constant_voltage(double const delta_t,
                                                     double const voltage)
{
  U_C -= (voltage - U_C) * std::expm1(-delta_t / (R * C));
  U = voltage;
  I = (U - U_C) / R;
}

void SeriesRC::evolve_one_time_step_linear_voltage(double const delta_t,
                                                   double const voltage)
{
  U_C -= (U - U_C) * std::expm1(-delta_t / (R * C));
  U_C += (voltage - U) / delta_t *
         (delta_t + R * C * std::expm1(-delta_t / (R * C)));
  U = voltage;
  I = (U - U_C) / R;
}

void SeriesRC::evolve_one_time_step_linear_power(double const delta_t,
                                                 double const power)
{
  std::ignore = delta_t;
  std::ignore = power;

  throw std::runtime_error("This function is not implemented.");
}

void SeriesRC::evolve_one_time_step_linear_load(double const delta_t,
                                                double const load)
{
  std::ignore = delta_t;
  std::ignore = load;

  throw std::runtime_error("This function is not implemented.");
}

std::size_t SeriesRC::evolve_one_time_step_constant_power(
    double const delta_t, double const power, std::string const &method)
{
  // TODO: if P is zero do constant current 0
  double const P          = power;
  double const ATOL       = 1.0e-14;
  double const RTOL       = 1.0e-14;
  std::size_t const MAXIT = 30;
  double const TOL        = std::abs(P) * RTOL + ATOL;
  size_t k = 0;
  while (true)
  {
    ++k;
    I = P / U;
    if (method.compare("FIXED_POINT") == 0)
      U = (R + delta_t / C) * I + U_C;
    else if (method.compare("NEWTON") == 0)
      U += ((R + delta_t / C) * P / U - U + U_C) /
           ((R + delta_t / C) * P / (U * U) + 1.0);
    else
      throw std::runtime_error("invalid method " + method);
    if (std::abs(P - U * I) < TOL)
      break;
    if (k >= MAXIT)
      throw std::runtime_error(method + " fail to converge within " +
                               std::to_string(MAXIT) + " iterations");
  }
  U_C += I * delta_t / C;
  return k;
}

ParallelRC::ParallelRC(boost::mpi::communicator const &comm,
                       boost::property_tree::ptree const &ptree)
    : EnergyStorageDevice(comm)
{
  R_series   = ptree.get<double>("series_resistance");
  R_parallel = ptree.get<double>("parallel_resistance");
  C          = ptree.get<double>("capacitance");
  U = ptree.get<double>("initial_voltage", 0.0);
  this->reset_voltage(U);
}

void ParallelRC::print_data(std::ostream &os) const
{
  os << boost::format("  %10.7f  %10.7f  %10.7f  \n") % I % U % U_C;
}

void ParallelRC::reset_voltage(double const voltage)
{
  U   = voltage;
  U_C = R_parallel / (R_series + R_parallel) * U;
  I   = U / (R_series + R_parallel);
}

void ParallelRC::reset_current(double const current)
{
  I = current;
  U = R_series * I + U_C;
}

void ParallelRC::reset(double const capacitor_voltage)
{
  U_C = capacitor_voltage;
  U   = R_series * I + U_C;
}

void ParallelRC::evolve_one_time_step_constant_current(double const delta_t,
                                                       double const current)
{
  U_C = R_parallel * current +
        (U_C - R_parallel * current) * std::exp(-delta_t / (R_parallel * C));
  I = current;
  U = R_series * I + U_C;
}

void ParallelRC::evolve_one_time_step_linear_current(double const delta_t,
                                                     double const current)
{
  U_C = R_parallel * I +
        (U_C - R_parallel * I) * std::exp(-delta_t / (R_parallel * C));
  U_C += R_parallel * (current - I) / delta_t *
         (delta_t + (R_parallel * C) * std::expm1(-delta_t / (R_parallel * C)));
  I = current;
  U = R_series * I + U_C;
}

void ParallelRC::evolve_one_time_step_constant_voltage(double const delta_t,
                                                       double const voltage)
{
  U_C -= (voltage * R_parallel / (R_series + R_parallel) - U_C) *
         std::expm1(-delta_t * (R_series + R_parallel) /
                    (R_series * R_parallel * C));
  U = voltage;
  I = (U - U_C) / R_series;
}

void ParallelRC::evolve_one_time_step_linear_voltage(double const delta_t,
                                                     double const voltage)
{
  U_C -= (U * R_parallel / (R_series + R_parallel) - U_C) *
         std::expm1(-delta_t * (R_series + R_parallel) /
                    (R_series * R_parallel * C));
  U_C += (voltage - U) / delta_t * R_parallel / (R_series + R_parallel) *
         (delta_t +
          (R_series * R_parallel * C) / (R_series + R_parallel) *
              std::expm1(-delta_t * (R_series + R_parallel) /
                         (R_series * R_parallel * C)));
  U = voltage;
  I = (U - U_C) / R_series;
}

void ParallelRC::evolve_one_time_step_linear_power(double const delta_t,
                                                   double const power)
{
  std::ignore = delta_t;
  std::ignore = power;

  throw std::runtime_error("This function is not implemented.");
}

void ParallelRC::evolve_one_time_step_linear_load(double const delta_t,
                                                  double const load)
{
  std::ignore = delta_t;
  std::ignore = load;

  throw std::runtime_error("This function is not implemented.");
}

void ParallelRC::evolve_one_time_step_constant_load(double const delta_t,
                                                    double const load)
{
  U_C *= std::exp(-delta_t * (1.0 + (R_series + load) / R_parallel) /
                  ((R_series + load) * C));
  I = -U_C / (R_series + load);
  U = U_C + R_series * I;
}

void ParallelRC::evolve_one_time_step_constant_power(double const delta_t,
                                                     double const power)
{
  evolve_one_time_step_constant_power(delta_t, power, "NEWTON");
}

std::size_t ParallelRC::evolve_one_time_step_constant_power(
    double const delta_t, double const power, std::string const &method)
{
  // TODO: if P is zero do constant current 0
  double const P          = power;
  double const ATOL       = 1.0e-14;
  double const RTOL       = 1.0e-14;
  std::size_t const MAXIT = 30;
  double const TOL        = std::abs(P) * RTOL + ATOL;
  size_t k = 0;
  while (true)
  {
    ++k;
    I = P / U;
    if (method.compare("FIXED_POINT") == 0)
      U = (R_series + R_parallel) * I +
          (U_C - R_parallel * I) * std::exp(-delta_t / (R_parallel * C));
    else if (method.compare("NEWTON") == 0)
      U += ((R_series +
             R_parallel * (1.0 - std::exp(-delta_t / (R_parallel * C)))) *
                P / U -
            U + U_C * std::exp(-delta_t / (R_parallel * C))) /
           ((R_series +
             R_parallel * (1.0 - std::exp(-delta_t / (R_parallel * C)))) *
                P / (U * U) +
            1.0);
    else
      throw std::runtime_error("invalid method " + method);
    if (std::abs(P - U * I) < TOL)
      break;
    if (k >= MAXIT)
      throw std::runtime_error(method + " fail to converge within " +
                               std::to_string(MAXIT) + " iterations");
  }
  U_C = U - R_series * I;
  return k;
}

} // end namespace
