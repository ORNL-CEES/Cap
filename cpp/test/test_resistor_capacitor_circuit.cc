/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#define BOOST_TEST_MODULE ResistorCapacitorCircuit
#define BOOST_TEST_MAIN
#include <cap/resistor_capacitor.h>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <string>
#include <tuple>
#include <cmath>
#include <iostream>


// This file contains the following tests:
//  - Series RC constant voltage
//  - Series RC constant current
//  - Series RC constant power
//  - Series RC constant load
//  - Parallel RC constant current
//  - Parallel RC constant voltage
//  - Parallel RC constant power
//  - Parallel RC constant load


double const R_SERIES   = 55.0e-3;
double const R_PARALLEL =  2.5e6;
double const C          =  3.0;
double const TOLERANCE  =  1.0e-8; // in percentage units
double const I          =  0.006;
double const U          =  2.1;
double const P          =  0.0017;

boost::property_tree::ptree initialize_database()
{
    boost::property_tree::ptree database;
    database.put("series_resistance"  , R_SERIES  );
    database.put("parallel_resistance", R_PARALLEL);
    database.put("capacitance"        , C         );
    return database;
}

void set_voltage(cap::SeriesRC &rc, double voltage)
{
  rc.U_C = voltage;
  rc.U   = rc.U_C;
  rc.I   = 0.0;
}

void set_current(cap::SeriesRC &rc, double current)
{
  rc.I = current;
  rc.U = rc.U_C + rc.R * rc.I;
}

void set_voltage(cap::ParallelRC &rc, double voltage)
{
  rc.U   = voltage;
  rc.U_C = rc.R_parallel / (rc.R_series + rc.R_parallel) * rc.U;
  rc.I   = rc.U / (rc.R_series + rc.R_parallel);
}

void set_current(cap::ParallelRC &rc, double current)
{
  rc.I = current;
  rc.U = rc.R_series * rc.I + rc.U_C;
}



BOOST_AUTO_TEST_CASE( test_series_rc_constant_voltage )
{
    double const TAU     = R_SERIES * C;
    double const DELTA_T = 0.1 * TAU;
    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc(boost::mpi::communicator(), initialize_database());

    // CHARGE
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, U * (1.0 - std::exp(-t/TAU)), TOLERANCE);
        rc.evolve_one_time_step_constant_voltage(DELTA_T, U);
    }

    // DISCHARGE
    set_voltage(rc, U);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, U * std::exp(-t/TAU), TOLERANCE);
        rc.evolve_one_time_step_constant_voltage(DELTA_T, 0.0);
    }
}



BOOST_AUTO_TEST_CASE( test_series_rc_constant_current )
{
    double const TAU     = R_SERIES * C;
    double const DELTA_T = 0.1 * TAU;

    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc(boost::mpi::communicator(), initialize_database());

    // CHARGE
    set_voltage(rc, 0.0);
    set_current(rc, I);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U, I * (R_SERIES + t / C), TOLERANCE);
        rc.evolve_one_time_step_constant_current(DELTA_T, I);
    }

    // DISCHARGE
    set_voltage(rc, U);
    set_current(rc, -I);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U, U - I * (R_SERIES + t / C), TOLERANCE);
        rc.evolve_one_time_step_constant_current(DELTA_T, -I);
    }

}



BOOST_AUTO_TEST_CASE( test_series_rc_constant_power )
{
    double const TAU     = R_SERIES * C;
    double const DELTA_T = 0.1 * TAU;

    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc_newton     (boost::mpi::communicator(), initialize_database());
    cap::SeriesRC rc_fixed_point(boost::mpi::communicator(), initialize_database());

    // CHARGE
    set_current(rc_newton, 0.0);
    set_voltage(rc_newton, U);
    set_current(rc_fixed_point, 0.0);
    set_voltage(rc_fixed_point, U);

    BOOST_FOREACH(double const & t, time)
    {
        std::ignore = t; 
        BOOST_CHECK_CLOSE(rc_newton.U,   rc_fixed_point.U,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.I,   rc_fixed_point.I,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.U_C, rc_fixed_point.U_C, TOLERANCE);
        rc_newton     .evolve_one_time_step_constant_power(DELTA_T, P, "NEWTON"     );
        rc_fixed_point.evolve_one_time_step_constant_power(DELTA_T, P, "FIXED_POINT");
    }

   BOOST_CHECK_THROW(rc_newton.evolve_one_time_step_constant_power(DELTA_T, P, "INVALID_ROOT_FINDING_METHOD"), std::runtime_error);

    // DISCHARGE
    set_current(rc_newton, 0.0);
    set_voltage(rc_newton, U);
    set_current(rc_fixed_point, 0.0);
    set_voltage(rc_fixed_point, U);
    BOOST_FOREACH(double const & t, time)
    {
        std::ignore = t; 
        BOOST_CHECK_CLOSE(rc_newton.U,   rc_fixed_point.U,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.I,   rc_fixed_point.I,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.U_C, rc_fixed_point.U_C, TOLERANCE);
        rc_newton     .evolve_one_time_step_constant_power(DELTA_T, -P, "NEWTON"     );
        rc_fixed_point.evolve_one_time_step_constant_power(DELTA_T, -P, "FIXED_POINT");
    }
}



BOOST_AUTO_TEST_CASE( test_series_rc_constant_load )
{
    double const TAU     = R_SERIES * C;
    double const DELTA_T = 0.1 * TAU;
    double const R_LOAD  = 5.0 * R_SERIES;

    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc(boost::mpi::communicator(), initialize_database());

    // DISCHARGE
    set_voltage(rc, U);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U, U * std::exp(- t / ((R_SERIES + R_LOAD) * C)) * (1.0 - ((t > 0.0) ? R_SERIES / (R_SERIES + R_LOAD) : 0.0)), TOLERANCE);
        rc.evolve_one_time_step_constant_load(DELTA_T, R_LOAD);
    }

}



BOOST_AUTO_TEST_CASE( test_parallel_rc_constant_current )
{
    double const TAU     = R_PARALLEL * C;
    double const DELTA_T = 0.1 * TAU;
    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += 0.1*TAU)
        time.push_back(t);

    cap::ParallelRC rc(boost::mpi::communicator(), initialize_database());
    rc.R_series = 0.0;

    // CHARGE
    set_voltage(rc, 0.0);
    set_current(rc, I);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, R_PARALLEL * I * (1.0 - std::exp(-t/(R_PARALLEL*C))), TOLERANCE);
        rc.evolve_one_time_step_constant_current(DELTA_T, I);
    }

    // RELAXATION
    set_voltage(rc,U  );
    set_current(rc,0.0);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, U * std::exp(-t/(R_PARALLEL*C)), TOLERANCE);
        rc.evolve_one_time_step_constant_current(DELTA_T, 0.0);
    }
}



BOOST_AUTO_TEST_CASE( test_parallel_rc_constant_voltage )
{
    double const TAU     = R_SERIES * C;
    double const DELTA_T = 0.1 * TAU;

    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::ParallelRC rc(boost::mpi::communicator(), initialize_database());

    // CHARGE
    set_voltage(rc, 0.0);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, 
            U * R_PARALLEL/(R_SERIES+R_PARALLEL)
                * (1.0 - std::exp(-t / ((R_SERIES*R_PARALLEL)/(R_SERIES+R_PARALLEL) * C))),
            TOLERANCE);
        rc.evolve_one_time_step_constant_voltage(DELTA_T, U);
    }

}



BOOST_AUTO_TEST_CASE( test_parallel_rc_constant_power )
{
    double const TAU     = R_SERIES * C;
    double const DELTA_T = 0.1 * TAU;

    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc_newton     (boost::mpi::communicator(), initialize_database());
    cap::SeriesRC rc_fixed_point(boost::mpi::communicator(), initialize_database());

    // CHARGE
    set_current(rc_newton, 0.0);
    set_voltage(rc_newton, U  );
    set_current(rc_fixed_point, 0.0);
    set_voltage(rc_fixed_point, U  );

    BOOST_FOREACH(double const & t, time)
    {
        std::ignore = t;
        BOOST_CHECK_CLOSE(rc_newton.U,   rc_fixed_point.U,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.I,   rc_fixed_point.I,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.U_C, rc_fixed_point.U_C, TOLERANCE);
        rc_newton     .evolve_one_time_step_constant_power(DELTA_T, P, "NEWTON"     );
        rc_fixed_point.evolve_one_time_step_constant_power(DELTA_T, P, "FIXED_POINT");
    }

    BOOST_CHECK_THROW(rc_newton.evolve_one_time_step_constant_power(DELTA_T, P, "INVALID_ROOT_FINDING_METHOD"), std::runtime_error);

    // DISCHARGE
    set_current(rc_newton, 0.0);
    set_voltage(rc_newton, U  );
    set_current(rc_fixed_point, 0.0);
    set_voltage(rc_fixed_point, U  );

    BOOST_FOREACH(double const & t, time)
    {
        std::ignore = t;
        BOOST_CHECK_CLOSE(rc_newton.U,   rc_fixed_point.U,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.I,   rc_fixed_point.I,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.U_C, rc_fixed_point.U_C, TOLERANCE);
        rc_newton     .evolve_one_time_step_constant_power(DELTA_T, -P, "NEWTON"     );
        rc_fixed_point.evolve_one_time_step_constant_power(DELTA_T, -P, "FIXED_POINT");
    }
}



BOOST_AUTO_TEST_CASE( test_parallel_rc_constant_load )
{
    double const TAU     = R_PARALLEL * C;
    double const DELTA_T = 0.1 * TAU;
    double const R_LOAD  = 5.0 * R_SERIES;
    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += 0.1*TAU)
        time.push_back(t);

    cap::ParallelRC rc(boost::mpi::communicator(), initialize_database());
    rc.R_series = 0.0;

    // DISCHARGE
    set_voltage(rc, U);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U, U * std::exp(- t * (1.0 + (R_SERIES + R_LOAD)/R_PARALLEL) / ((R_SERIES + R_LOAD) * C)) * (1.0 - ((t > 0.0) ? R_SERIES * (1.0 + (R_SERIES + R_LOAD)/R_PARALLEL) / (R_SERIES + R_LOAD) : 0.0)), TOLERANCE);
        rc.evolve_one_time_step_constant_load(DELTA_T, R_LOAD);
    }
}
