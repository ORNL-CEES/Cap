/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#define BOOST_TEST_MODULE ExactTransientSolution
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/mp_values.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

namespace cap {

void verification_problem(std::shared_ptr<cap::EnergyStorageDevice> dev)
{
    // gold vs computed
    double const charge_current = 5e-3;
    double const charge_time    = 0.01;
    double const time_step      = 1e-4;
    double const epsilon        = time_step * 1.0e-4;

    unsigned int pos = 0;
    double computed_voltage;
    double const percent_tolerance = 1e-6;
    std::vector<double> gold_solution(10);
    gold_solution[0] = 1.725914356067658e-01;
    gold_solution[1] = 1.802025636145941e-01;
    gold_solution[2] = 1.859326352495181e-01;
    gold_solution[3] = 1.905978440188036e-01;
    gold_solution[4] = 1.946022119085378e-01;
    gold_solution[5] = 1.981601232287249e-01;
    gold_solution[6] = 2.013936650249285e-01;
    gold_solution[7] = 2.043807296399895e-01;
    gold_solution[8] = 2.071701713934283e-01;
    gold_solution[9] = 2.097979282542038e-01;
    for (double time = 0.0; time <= charge_time+epsilon; time += time_step)
    {
        dev->evolve_one_time_step_constant_current(time_step, charge_current);
        dev->get_voltage(computed_voltage);
        if ((std::abs(time+time_step-1e-3) < 1e-7) || (std::abs(time+time_step-2e-3) < 1e-7) ||
            (std::abs(time+time_step-3e-3) < 1e-7) || (std::abs(time+time_step-4e-3) < 1e-7) ||
            (std::abs(time+time_step-5e-3) < 1e-7) || (std::abs(time+time_step-6e-3) < 1e-7) ||
            (std::abs(time+time_step-7e-3) < 1e-7) || (std::abs(time+time_step-8e-3) < 1e-7) ||
            (std::abs(time+time_step-9e-3) < 1e-7) || (std::abs(time+time_step-10e-3) < 1e-7))
        {
          BOOST_CHECK_CLOSE(computed_voltage, gold_solution[pos], percent_tolerance);
          ++pos;
        }
    }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_exact_transient_solution )
{
    // parse input file
    boost::property_tree::ptree device_database;
    boost::property_tree::info_parser::read_info("super_capacitor.info", device_database);

    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::EnergyStorageDevice::build(boost::mpi::communicator(), device_database);

    cap::verification_problem(device);
}    
