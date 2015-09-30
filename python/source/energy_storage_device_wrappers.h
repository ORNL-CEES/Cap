#include <cap/energy_storage_device.h>
#include <cap/equivalent_circuit.h>
#include <boost/python/wrapper.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <map>

namespace pycap {

struct EnergyStorageDeviceWrap : cap::EnergyStorageDevice, boost::python::wrapper<cap::EnergyStorageDevice>
{
    void evolve_one_time_step_constant_current(double const time_step, double const constant_current);
    void evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage);
    void evolve_one_time_step_constant_power  (double const time_step, double const constant_power  );
    void evolve_one_time_step_constant_load   (double const time_step, double const constant_load   );
    void evolve_one_time_step_changing_current(double const time_step, double const changing_current);
    void evolve_one_time_step_changing_voltage(double const time_step, double const changing_voltage);
    void evolve_one_time_step_changing_power  (double const time_step, double const changing_power  );
    void evolve_one_time_step_changing_load   (double const time_step, double const changing_load   );
};

double get_current(cap::EnergyStorageDevice const & device);
double get_voltage(cap::EnergyStorageDevice const & device);

std::shared_ptr<cap::EnergyStorageDevice> build_energy_storage_device(boost::python::object & python_object);

// DEPRECATED
struct ElectrochemicalImpedanceSpectroscopyData {
    void impedance_spectroscopy(boost::python::object & python_device, boost::python::object & python_database);
    void measure_impedance(boost::python::object & python_device, boost::python::object & python_database);
    boost::python::list get_frequency() const;
    boost::python::list get_complex_impedance() const;
    void clear();
    std::map<double,std::complex<double>> data;
};

std::shared_ptr<boost::property_tree::ptree> compute_equivalent_circuit(boost::python::object & python_object);

} // end namespace pycap

