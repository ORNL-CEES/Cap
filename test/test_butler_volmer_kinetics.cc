#define BOOST_TEST_MODULE ButlerVolmerKinetics
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <cmath>
#include <iostream>


namespace cap {
class SeriesRC
{
public:
    SeriesRC(double const resistance, double const capacitance, double const initial_capacitor_voltage = 0.0)
        : R(resistance)
        , C(capacitance)
        , U_C(initial_capacitor_voltage) 
        , U(U_C)
        , I(0.0) 
    { }
    void reset(double const capacitor_voltage)
    {
        U_C = capacitor_voltage;
        U = U_C + R * I;
    }
    void evolve_one_time_step_constant_current(double const delta_t, double const constant_current)
    {
        I = constant_current; 
        U_C += I * delta_t / C;
        U = R * I + U_C;
    }
    void evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage)
    {
        U = constant_voltage; 
        U_C = U - (U - U_C) * std::exp(-delta_t/(R*C));
        I = U - U_C / R;
    }
    std::size_t evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method = "NEWTON")
    {
        // TODO: if P is zero do constant current 0
        double const P = constant_power;
        double const      ATOL  = 1.0e-14;
        double const      RTOL  = 1.0e-14;
        std::size_t const MAXIT = 30;
        double const TOL = std::abs(P) * RTOL + ATOL;
        size_t k = 0;
        while (true)
        {
            ++k;
            I = P / U;
            if (method.compare("FIXED_POINT") == 0)
                U = (R + delta_t / C) * I + U_C;
            else if (method.compare("NEWTON") == 0)
                U += ((R + delta_t / C) * P / U - U + U_C) / ((R + delta_t / C) * P / (U * U) + 1.0);
            else
                throw std::runtime_error("invalid method "+method);
            if (std::abs(P - U * I) < TOL)
                break;
            if (k >= MAXIT)
                throw std::runtime_error(method+" fail to converge within "+std::to_string(MAXIT)+" iterations");
        }
        U_C += I * delta_t / C;
        return k;
    }
    
    double R;
    double C;
    double U_C;
    double U;
    double I;

};
} // end namespace


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
