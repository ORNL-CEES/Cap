#include <cap/resistor_capacitor.h>
#include <cmath>
#include <stdexcept>

namespace cap {

SeriesRC::
SeriesRC(double const resistance, double const capacitance, double const initial_capacitor_voltage)
    : R(resistance)
    , C(capacitance)
    , U_C(initial_capacitor_voltage) 
    , U(U_C)
    , I(0.0) 
{ }



void
SeriesRC::
reset(double const capacitor_voltage)
{
    U_C = capacitor_voltage;
    U = U_C + R * I;
}



void
SeriesRC::
evolve_one_time_step_constant_current(double const delta_t, double const constant_current)
{
    I = constant_current; 
    U_C += I * delta_t / C;
    U = R * I + U_C;
}



void
SeriesRC::
evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage)
{
    U = constant_voltage; 
    U_C = U - (U - U_C) * std::exp(-delta_t/(R*C));
    I = U - U_C / R;
}



std::size_t
SeriesRC::
evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method)
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

} // end namespace
