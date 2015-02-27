#ifndef CAP_RESISTOR_CAPACITOR_H
#define CAP_RESISTOR_CAPACITOR_H

#include <string>

namespace cap {

class SeriesRC
{
public:
    SeriesRC(double const resistance, double const capacitance, double const initial_capacitor_voltage = 0.0);
    void reset(double const capacitor_voltage);
    void evolve_one_time_step_constant_current(double const delta_t, double const constant_current);
    void evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage);
    std::size_t evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method = "NEWTON");
    
    double R;
    double C;
    double U_C;
    double U;
    double I;
};

} // end namespace

#endif // CAP_RESISTOR_CAPACITOR_H
