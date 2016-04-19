/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE TestImpedanceSpectroscopy
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <cap/utils.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_fft_real.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <iterator>
#include <algorithm>

namespace cap
{

std::function<std::vector<std::tuple<double, std::complex<double>>>(
    std::shared_ptr<boost::property_tree::ptree const>)>
get_compute_exact(
    std::shared_ptr<boost::property_tree::ptree const> dev_database,
    std::shared_ptr<boost::property_tree::ptree const> eis_database)
{
  std::string const device_type = dev_database->get<std::string>("type");
  double const series_resistance =
      dev_database->get<double>("series_resistance");
  double const parallel_resistance =
      dev_database->get<double>("parallel_resistance");
  double const capacitance = dev_database->get<double>("capacitance");

  std::vector<int> const harmonics =
      cap::to_vector<int>(eis_database->get<std::string>("harmonics"));
  double const pi = boost::math::constants::pi<double>();

  if (device_type.compare("SeriesRC") == 0)
  {
    return [series_resistance, capacitance, pi, harmonics](
        std::shared_ptr<boost::property_tree::ptree const> database)
    {
      std::vector<std::tuple<double, std::complex<double>>> results;
      double const frequency = database->get<double>("frequency");
      for (int k : harmonics)
        results.emplace_back(std::make_tuple(
            k * frequency,
            series_resistance +
                1.0 / std::complex<double>(0.0, capacitance * 2.0 * pi * k *
                                                    frequency)));
      return results;
    };
  }
  else if (device_type.compare("ParallelRC") == 0)
  {
    return [series_resistance, parallel_resistance, capacitance, pi, harmonics](
        std::shared_ptr<boost::property_tree::ptree const> database)
    {
      std::vector<std::tuple<double, std::complex<double>>> results;
      double const frequency = database->get<double>("frequency");
      for (int k : harmonics)
        results.emplace_back(std::make_tuple(
            k * frequency, series_resistance +
                               parallel_resistance /
                                   std::complex<double>(
                                       1.0, parallel_resistance * capacitance *
                                                2.0 * pi * k * frequency)));
      return results;
    };
  }
  else
  {
    throw std::runtime_error("invalid device type " + device_type);
  }
}

void scan(std::shared_ptr<cap::EnergyStorageDevice> dev,
          std::shared_ptr<boost::property_tree::ptree const> dev_database,
          std::shared_ptr<boost::property_tree::ptree const> eis_database,
          std::ostream &os = std::cout)
{
  double const frequency_upper_limit =
      eis_database->get<double>("frequency_upper_limit");
  double const frequency_lower_limit =
      eis_database->get<double>("frequency_lower_limit");
  int const steps_per_decade = eis_database->get<int>("steps_per_decade");
  double const pi = boost::math::constants::pi<double>();
  std::vector<int> const harmonics =
      cap::to_vector<int>(eis_database->get<std::string>("harmonics"));
  double const percent_tolerance =
      eis_database->get<double>("percent_tolerance");

  auto compute_exact = get_compute_exact(dev_database, eis_database);

  double expe_frequency;
  std::complex<double> expe_impedance;
  double theo_frequency;
  std::complex<double> theo_impedance;

  std::shared_ptr<boost::property_tree::ptree> tmp =
      std::make_shared<boost::property_tree::ptree>(*eis_database);

  os << "# impedance Z(f) = R + i X \n";
  os << boost::format("# %22s  %22s  %22s  %22s  %22s  \n") %
            "frequency_f_[Hz]" % "resistance_R_[ohm]" % "reactance_X_[ohm]" %
            "magnitude_|Z|_[ohm]" % "phase_arg(Z)_[degree]";

  for (double frequency = frequency_upper_limit;
       frequency >= frequency_lower_limit;
       frequency /= std::pow(10.0, 1.0 / steps_per_decade))
  {
    tmp->put("frequency", frequency);
    auto expe_results = measure_impedance(dev, tmp);
    auto theo_results = compute_exact(tmp);
    for (std::size_t k = 0; k < harmonics.size(); ++k)
    {
      auto it = expe_results.begin();
      std::advance(it, k);
      std::tie(expe_frequency, expe_impedance) = *it;
      std::tie(theo_frequency, theo_impedance) = theo_results[k];

      BOOST_CHECK_CLOSE(expe_frequency, theo_frequency, percent_tolerance);
      BOOST_CHECK_CLOSE(expe_impedance.real(), theo_impedance.real(),
                        percent_tolerance);
      BOOST_CHECK_CLOSE(expe_impedance.imag(), theo_impedance.imag(),
                        percent_tolerance);
      BOOST_CHECK_CLOSE(std::abs(expe_impedance), std::abs(theo_impedance),
                        percent_tolerance);
      BOOST_CHECK_CLOSE(std::arg(expe_impedance), std::arg(theo_impedance),
                        percent_tolerance);

      os << boost::format("  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  "
                          "%22.15e  %22.15e  %22.15e  %22.15e  %22.15e \n") %
                expe_frequency % expe_impedance.real() % expe_impedance.imag() %
                std::abs(expe_impedance) %
                (std::arg(expe_impedance) * 180.0 / pi) % theo_frequency %
                theo_impedance.real() % theo_impedance.imag() %
                std::abs(theo_impedance) %
                (std::arg(theo_impedance) * 180.0 / pi);
    }
  }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE(test_impedance_spectroscopy)
{
  // parse input file
  std::shared_ptr<boost::property_tree::ptree> input_database =
      std::make_shared<boost::property_tree::ptree>();
  boost::property_tree::info_parser::read_info(
      "input_impedance_spectroscopy.info", *input_database);

  // build an energy storage system
  std::shared_ptr<boost::property_tree::ptree> device_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database->get_child("device"));
  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(boost::mpi::communicator(),
                                      *device_database);

  // measure its impedance
  std::fstream fout;
  fout.open("computed_vs_exact_impedance_spectroscopy_data", std::fstream::out);

  std::shared_ptr<boost::property_tree::ptree> impedance_spectroscopy_database =
      std::make_shared<boost::property_tree::ptree>(
          input_database->get_child("impedance_spectroscopy"));
  cap::scan(device, device_database, impedance_spectroscopy_database, fout);

  fout.close();
}
