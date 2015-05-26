#ifndef CAP_RESISTOR_CAPACITOR_H
#define CAP_RESISTOR_CAPACITOR_H

#include <cap/energy_storage_device.h>
#include <string>

namespace cap {

class SeriesRC : public EnergyStorageDevice
{
public:
    SeriesRC(std::shared_ptr<Parameters const> );
    void print_data(std::ostream & os) const override;
    void reset_voltage(double const voltage) override;
    void reset_current(double const current) override;
    void evolve_one_time_step_constant_current(double const delta_t, double const constant_current) override;
    void evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage) override;
    void evolve_one_time_step_constant_power  (double const delta_t, double const constant_power  ) override;
    void evolve_one_time_step_constant_load   (double const delta_t, double const constant_load   ) override;
    void evolve_one_time_step_changing_current(double const delta_t, double const changing_current) override;
    void evolve_one_time_step_changing_voltage(double const delta_t, double const changing_voltage) override;
    inline void get_voltage(double & voltage) const override { voltage = U; }
    inline void get_current(double & current) const override { current = I; }
    // TODO:
    std::size_t evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method = "NEWTON");
    void reset(double const capacitor_voltage);
    
    double R;
    double C;
    double U_C;
    double U;
    double I;
};

class ParallelRC : public EnergyStorageDevice
{
public:
    ParallelRC(std::shared_ptr<Parameters const> );
    void print_data(std::ostream & os) const override;
    void reset_voltage(double const voltage) override;
    void reset_current(double const current) override;
    void evolve_one_time_step_constant_current(double const delta_t, double const constant_current) override;
    void evolve_one_time_step_constant_voltage(double const delta_t, double const constant_voltage) override;
    void evolve_one_time_step_constant_power  (double const delta_t, double const constant_power  ) override;
    void evolve_one_time_step_constant_load   (double const delta_t, double const constant_load   ) override;
    void evolve_one_time_step_changing_current(double const delta_t, double const changing_current) override;
    void evolve_one_time_step_changing_voltage(double const delta_t, double const constant_voltage) override;
    inline void get_voltage(double & voltage) const override { voltage = U; }
    inline void get_current(double & current) const override { current = I; }
    // TODO:
    std::size_t evolve_one_time_step_constant_power(double const delta_t, double const constant_power, std::string const & method = "NEWTON");
    void reset(double const capacitor_voltage);

    double R_parallel;
    double C;
    double U_C;
    double U;
    double I;
    double R_series;
};


} // end namespace

#endif // CAP_RESISTOR_CAPACITOR_H
