/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/electrochemical_impedance_spectroscopy.h>
#include <cap/utils.h>
#include <boost/assert.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define REAL(z, i) ((z)[2 * (i)])
#define IMAG(z, i) ((z)[2 * (i) + 1])
#include <iostream>

namespace cap
{

std::vector<std::complex<double>> compute_fft(double const *input,
                                              std::size_t const n)
{
  std::vector<double> data(input, input + n);
  gsl_fft_real_radix2_transform(&(data[0]), 1, n);
  std::vector<double> unpacked_data(2 * n);
  gsl_fft_halfcomplex_radix2_unpack(&(data[0]), &(unpacked_data[0]), 1, n);
  std::vector<std::complex<double>> output(n / 2 + 1);
  for (std::size_t k = 0; k < n / 2 + 1; ++k)
    output[k] =
        std::complex<double>(REAL(unpacked_data, k), IMAG(unpacked_data, k));
  return output;
}

std::vector<double> compute_fft_frequencies(std::size_t const n,
                                            double const frequency)
{
  std::vector<double> fft_frequencies(n / 2 + 1);
  for (std::size_t k = 0; k < n / 2 + 1; ++k)
    fft_frequencies[k] = k * frequency;
  return fft_frequencies;
}

std::vector<double>
compute_signal_to_noise_ratio(std::vector<std::complex<double>> fft_data,
                              std::vector<int> excited_harmonics,
                              std::vector<int> unexcited_harmonics)
{
  BOOST_ASSERT(fft_data.size() ==
               excited_harmonics.size() + unexcited_harmonics.size() + 1);
  double dc_power = std::norm(fft_data[0]);
  double ac_power = 0.0;
  for (int k : excited_harmonics)
    ac_power += std::norm(fft_data[k]);
  double noise_power = 0.0;
  for (int k : unexcited_harmonics)
    noise_power += std::norm(fft_data[k]);
  double total_power = dc_power + ac_power + noise_power;
  std::ignore = total_power;
  // TODO:  definition is a bit shaky for multiple frequency
  // should compute for each individual excited freq vs nearby unex freq
  std::vector<double> ratios;
  for (std::size_t k : excited_harmonics)
  {
    ratios.emplace_back(std::norm(fft_data[k]) / noise_power);
  }
  return ratios;
}

std::function<void(std::shared_ptr<cap::EnergyStorageDevice>, double, double)>
get_evolve_one_time_step(boost::property_tree::ptree const &ptree)
{
  // clang-format off
  std::vector<int> const harmonics     = cap::to_vector<int>(ptree.get<std::string>("harmonics"));
  std::vector<double> const amplitudes = cap::to_vector<double>(ptree.get<std::string>("amplitudes"));
  std::vector<double> const phases     = cap::to_vector<double>(ptree.get<std::string>("phases"));
  // clang-format on
  BOOST_ASSERT(harmonics.size() == amplitudes.size());
  BOOST_ASSERT(harmonics.size() == phases.size());
  double const frequency = ptree.get<double>("frequency");
  auto compute_ac_excitation_signal = [frequency, harmonics, amplitudes,
                                       phases](double time)
  {
    double excitation_signal = 0.0;
    for (std::size_t k = 0; k < harmonics.size(); ++k)
      excitation_signal +=
          amplitudes[k] *
          boost::math::sin_pi(2 * harmonics[k] * frequency * time + phases[k]);
    return excitation_signal;
  };

  std::string const mode = ptree.get<std::string>("mode");
  if (mode.compare("galvanostatic") != 0)
  {
    return [compute_ac_excitation_signal](
        std::shared_ptr<cap::EnergyStorageDevice> device, double time,
        double time_step)
    {
      device->evolve_one_time_step_linear_current(
          time_step, compute_ac_excitation_signal(time));
    };
  }
  else if (mode.compare("potentiostatic") != 0)
  {
    double const dc_voltage = ptree.get<double>("dc_voltage");
    return [dc_voltage, compute_ac_excitation_signal](
        std::shared_ptr<cap::EnergyStorageDevice> device, double time,
        double time_step)
    {
      device->evolve_one_time_step_linear_voltage(
          time_step, dc_voltage + compute_ac_excitation_signal(time));
    };
  }
  else
  {
    throw std::runtime_error("invalide mode " + mode + " for EIS measurement");
  }
}

std::map<double, std::complex<double>>
measure_impedance(std::shared_ptr<cap::EnergyStorageDevice> device,
                  boost::property_tree::ptree const &ptree)
{
  std::vector<int> const harmonics =
      cap::to_vector<int>(ptree.get<std::string>("harmonics"));
  double const frequency = ptree.get<double>("frequency");
  int const cycles = ptree.get<int>("cycles");
  int const ignore_cycles = ptree.get<int>("ignore_cycles");
  int const steps_per_cycle = ptree.get<int>("steps_per_cycle");
  auto evolve_one_time_step = get_evolve_one_time_step(ptree);

  // apply excitation signal and measure response
  std::vector<double> time(cycles * steps_per_cycle);
  std::vector<double> current(cycles * steps_per_cycle);
  std::vector<double> voltage(cycles * steps_per_cycle);
  double const time_step = 1.0 / frequency / steps_per_cycle;
  for (int step = 0; step < cycles * steps_per_cycle; ++step)
  {
    time[step] = (step + 1) * time_step;
    evolve_one_time_step(device, time[step], time_step);
    device->get_current(current[step]);
    device->get_voltage(voltage[step]);
  }

  // compute discrete fourrier transforms
  std::vector<std::complex<double>> fft_current =
      compute_fft(&(current[ignore_cycles * steps_per_cycle]),
                  (cycles - ignore_cycles) * steps_per_cycle);
  std::vector<std::complex<double>> fft_voltage =
      compute_fft(&(voltage[ignore_cycles * steps_per_cycle]),
                  (cycles - ignore_cycles) * steps_per_cycle);
  std::vector<double> fft_frequencies =
      compute_fft_frequencies((cycles - ignore_cycles) * steps_per_cycle,
                              frequency / (cycles - ignore_cycles));

  //
  int const n = (cycles - ignore_cycles) * steps_per_cycle;
  BOOST_ASSERT(fft_current.size() == static_cast<std::size_t>(n / 2 + 1));
  BOOST_ASSERT(fft_voltage.size() == static_cast<std::size_t>(n / 2 + 1));
  BOOST_ASSERT(fft_frequencies.size() == static_cast<std::size_t>(n / 2 + 1));
  std::vector<int> excited_harmonics;
  std::vector<int> unexcited_harmonics;
  for (int i = 1; i <= n / 2; ++i)
    if ((std::find(harmonics.begin(), harmonics.end(),
                   i / (cycles - ignore_cycles)) != harmonics.end()) &&
        (i % (cycles - ignore_cycles) == 0))
      excited_harmonics.push_back(i);
    else
      unexcited_harmonics.push_back(i);

  // check signal-to-noise ratio
  std::vector<double> current_signal_to_noise_ratio =
      compute_signal_to_noise_ratio(fft_current, excited_harmonics,
                                    unexcited_harmonics);
  std::vector<double> voltage_signal_to_noise_ratio =
      compute_signal_to_noise_ratio(fft_voltage, excited_harmonics,
                                    unexcited_harmonics);
  std::cout << "current signal to noise ratio = ";
  for (auto const &ratio : current_signal_to_noise_ratio)
    std::cout << ratio << "  ";
  std::cout << "\n";
  std::cout << "voltage signal to noise ratio = ";
  for (auto const &ratio : voltage_signal_to_noise_ratio)
    std::cout << ratio << "  ";
  std::cout << "\n";
  // TODO: do we want to assert something or give a warning if ratio is not good
  // enough?

  // return complex impedance
  std::map<double, std::complex<double>> impedance;
  for (int const &k : excited_harmonics)
  {
    impedance.emplace(fft_frequencies[k], fft_voltage[k] / fft_current[k]);
  }
  return impedance;
}

std::map<double, std::complex<double>>
impedance_spectroscopy(std::shared_ptr<cap::EnergyStorageDevice> device,
                       boost::property_tree::ptree const &database)
{
  // clang-format off
  double const frequency_upper_limit = database.get<double>("frequency_upper_limit");
  double const frequency_lower_limit = database.get<double>("frequency_lower_limit");
  int const steps_per_decade         = database.get<int>("steps_per_decade");
  // clang-format on
  boost::property_tree::ptree tmp_database(database);
  std::map<double, std::complex<double>> data;
  for (double frequency = frequency_upper_limit;
       frequency >= frequency_lower_limit;
       frequency /= std::pow(10.0, 1.0 / steps_per_decade))
  {
    tmp_database.put("frequency", frequency);
    auto const impedance = measure_impedance(device, tmp_database);
    data.insert(impedance.begin(), impedance.end());
  }
  return data;
}

} // end namespace cap
