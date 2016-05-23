/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE TestPostProcessor

#include "main.cc"

#include <cap/post_processor.h>
#include <cap/utils.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

BOOST_AUTO_TEST_CASE(test_compute_energy)
{
  // We need to use a large tolerance because the error increases at each time
  // step.
  double const tolerance = 1e-6;

  double const pi = M_PI;
  double const initial_time = 0.0;
  double const final_time = 2.0 * pi;
  std::size_t const n = 10001;
  std::vector<double> time(n);
  std::vector<double> power(n);

  // INITIALIZE SOLUTION
  for (std::size_t i = 0; i < n; ++i)
  {
    time[i] = initial_time +
              static_cast<double>(i) / (n - 1) * (final_time - initial_time);
    power[i] = std::cos(time[i]);
  }

  // COMPUTE ENERGY
  std::vector<double> energy(n);
  cap::compute_energy(time, power, energy);

  // CHECK THE ANSWER
  std::vector<double> exact(n);
  std::vector<double> error(n);
  exact[0] = 0.0;
  for (std::size_t i = 0; i < n; ++i)
  {
    exact[i] = std::sin(time[i]);
    error[i] = energy[i] - exact[i];
    BOOST_REQUIRE(std::abs(error[i]) < tolerance);
  }
}
