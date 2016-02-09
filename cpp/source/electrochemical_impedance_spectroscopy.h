/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#ifndef CAP_ELECTROCHEMICAL_IMPEDANCE_SPECTROSCOPY_H
#define CAP_ELECTROCHEMICAL_IMPEDANCE_SPECTROSCOPY_H
#include <cap/energy_storage_device.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <vector>
#include <map>
#include <complex>

namespace cap
{

std::vector<std::complex<double>> compute_fft(double const *data,
                                              size_t const n);

std::vector<double> compute_fft_frequencies(size_t const n, double const f);

std::map<double, std::complex<double>>
measure_impedance(std::shared_ptr<cap::EnergyStorageDevice> dev,
                  std::shared_ptr<boost::property_tree::ptree const> database);

std::map<double, std::complex<double>> impedance_spectroscopy(
    std::shared_ptr<cap::EnergyStorageDevice> dev,
    std::shared_ptr<boost::property_tree::ptree const> database);

} // end namespace cap
#endif
