#define BOOST_TEST_MODULE ButlerVolmerKinetics
#define BOOST_TEST_MAIN
#include <cap/resistor_capacitor.h>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <cmath>
#include <iostream>



double const I         =  0.006;
double const U         =  2.1;
double const R         = 55.0;
double const C         =  3.0;
double const P         =  0.0017;
double const TAU       = R * C;
double const DELTA_T   = 0.1 * TAU;
double const TOLERANCE = 1.0e-8;    // in percentage units

BOOST_AUTO_TEST_CASE( test_series_rc_constant_voltage )
{
    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc(R, C);

    // CHARGE
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, U * (1.0 - std::exp(-t/TAU)), TOLERANCE);
        rc.evolve_one_time_step_constant_voltage(DELTA_T, U);
    }

    // DISCHARGE
    rc.reset(U);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U_C, U * std::exp(-t/TAU), TOLERANCE);
        rc.evolve_one_time_step_constant_voltage(DELTA_T, 0.0);
    }
}

BOOST_AUTO_TEST_CASE( test_series_rc_constant_current )
{
    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc(R, C);

    // CHARGE
    rc.I = I;
    rc.reset(0.0);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U, I * (R + t / C), TOLERANCE);
        rc.evolve_one_time_step_constant_current(DELTA_T, I);
    }

    // DISCHARGE
    rc.I = -I;
    rc.reset(U);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc.U, U - I * (R + t / C), TOLERANCE);
        rc.evolve_one_time_step_constant_current(DELTA_T, -I);
    }

}

BOOST_AUTO_TEST_CASE( test_series_rc_constant_power )
{
    std::vector<double> time;
    for (double t = 0.0; t <= 5.0 * TAU; t += DELTA_T)
        time.push_back(t);

    cap::SeriesRC rc_newton     (R, C);
    cap::SeriesRC rc_fixed_point(R, C);

    // CHARGE
    rc_newton     .reset(U);
    rc_fixed_point.reset(U);

    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc_newton.U,   rc_fixed_point.U,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.I,   rc_fixed_point.I,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.U_C, rc_fixed_point.U_C, TOLERANCE);
        rc_newton     .evolve_one_time_step_constant_power(DELTA_T, P, "NEWTON"     );
        rc_fixed_point.evolve_one_time_step_constant_power(DELTA_T, P, "FIXED_POINT");
    }

   BOOST_CHECK_THROW(rc_newton.evolve_one_time_step_constant_power(DELTA_T, P, "INVALID_ROOT_FINDING_METHOD"), std::runtime_error);

    // DISCHARGE
    rc_newton     .reset(U);
    rc_fixed_point.reset(U);
    BOOST_FOREACH(double const & t, time)
    {
        BOOST_CHECK_CLOSE(rc_newton.U,   rc_fixed_point.U,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.I,   rc_fixed_point.I,   TOLERANCE);
        BOOST_CHECK_CLOSE(rc_newton.U_C, rc_fixed_point.U_C, TOLERANCE);
        rc_newton     .evolve_one_time_step_constant_power(DELTA_T, -P, "NEWTON"     );
        rc_fixed_point.evolve_one_time_step_constant_power(DELTA_T, -P, "FIXED_POINT");
    }
}
