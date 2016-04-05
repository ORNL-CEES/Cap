/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/post_processor.templates.h>
#include <boost/numeric/odeint.hpp>

namespace cap
{

template class PostprocessorParameters<2>;
template class Postprocessor<2>;
template class PostprocessorParameters<3>;
template class Postprocessor<3>;

template class SuperCapacitorPostprocessorParameters<2>;
template class SuperCapacitorPostprocessor<2>;
template class SuperCapacitorPostprocessorParameters<3>;
template class SuperCapacitorPostprocessor<3>;

namespace internal
{
class energy_odeint
{
public:
  energy_odeint(std::vector<double> const &time,
                std::vector<double> const &power);

  void operator()(std::vector<double> const &, std::vector<double> &dxdt,
                  double const current_time);

private:
  std::vector<double> const time;
  std::vector<double> const power;
};

energy_odeint::energy_odeint(std::vector<double> const &time,
                             std::vector<double> const &power)
    : time(time), power(power)
{
}

void energy_odeint::operator()(std::vector<double> const &,
                               std::vector<double> &dxdt,
                               double const current_time)
{
  unsigned int pos = 0;
  for (unsigned int i = 0; i < time.size(); ++i, ++pos)
  {
    if (time[i] >= current_time)
      break;
  }
  if (pos == 0)
    dxdt[0] = power[0];
  else
    dxdt[0] = power[pos - 1] +
              (current_time - time[pos - 1]) / (time[pos] - time[pos - 1]) *
                  (power[pos] - power[pos - 1]);
}

void evaluate_energy(std::vector<double> const &computed_time,
                     std::vector<double> const &computed_energy,
                     std::vector<double> const &time,
                     std::vector<double> &energy)
{
  energy.clear();

  for (auto current_time : time)
  {
    unsigned int pos = 0;
    for (unsigned int i = 0; i < computed_time.size(); ++i, ++pos)
    {
      if (computed_time[i] >= current_time)
        break;
    }
    if (pos == 0)
      energy.push_back(computed_energy[0]);
    else
      energy.push_back(computed_energy[pos - 1] +
                       (current_time - computed_time[pos - 1]) /
                           (computed_time[pos] - computed_time[pos - 1]) *
                           (computed_energy[pos] - computed_energy[pos - 1]));
  }
}
}

void compute_energy(std::vector<double> const &time,
                    std::vector<double> const &power,
                    std::vector<double> &energy)

{
  bool const valid_input =
      (time.size() == power.size()) && (time.size() == energy.size());
  if (!valid_input)
    throw std::runtime_error("invalid input");
  internal::energy_odeint ode(time, power);
  std::vector<double> computed_energy;
  std::vector<double> computed_time;
  std::vector<double> x0(1);
  x0[0] = 0.;
  // Find the smallest time step
  double dt = time[1] - time[0];
  if (time.size() > 2)
    for (unsigned int i = 2; i < time.size(); ++i)
      if (dt > time[i] - time[i - 1])
        dt = time[i] - time[i - 1];

  double const abs_tol = 1e-9;
  double const rel_tol = 1e-9;
  boost::numeric::odeint::integrate_const(
      boost::numeric::odeint::make_dense_output(
          abs_tol, rel_tol,
          boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>>()),
      ode, x0, time[0], time.back(), dt / 2.,
      [&](std::vector<double> const &x, double const &t)
      {
        computed_energy.push_back(x[0]);
        computed_time.push_back(t);
      });

  internal::evaluate_energy(computed_time, computed_energy, time, energy);
}

void extract_duration_and_average_power(
    std::vector<std::string> const &capacitor_state,
    std::vector<double> const &time, std::vector<double> const &energy,
    std::vector<double> &duration, std::vector<double> &average_power)
{
  bool const valid_input = (time.size() == energy.size()) && duration.empty() &&
                           average_power.empty();
  if (!valid_input)
    throw std::runtime_error("invalid input");
  std::vector<std::string>::const_iterator it     = capacitor_state.begin();
  std::vector<std::string>::const_iterator end_it = capacitor_state.end();
  std::size_t first = 0;
  while (it != end_it)
  {
    auto same = [&it](std::string const &o)
    {
      return it->compare(o) == 0;
    };
    std::vector<std::string>::const_iterator next =
        std::find_if_not(it, end_it, same);
    std::size_t last = first + std::distance(it, next);
    duration.push_back(time[last - 1] - time[first]);
    average_power.push_back((
        first != last - 1 ? (energy[last - 1] - energy[first]) / duration.back()
                          : 0.0));
    it    = next;
    first = last;
  }
}

} // end namespace cap
