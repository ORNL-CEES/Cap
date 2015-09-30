#include <pycap/energy_storage_device_wrappers.h>
#ifdef WITH_GSL
#include <cap/electrochemical_impedance_spectroscopy.h>
#endif

namespace pycap {

void EnergyStorageDeviceWrap::evolve_one_time_step_constant_current(double const time_step, double const constant_current)
{
    this->get_override("evolve_one_time_step_constant_current")(time_step, constant_current);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage)
{
    this->get_override("evolve_one_time_step_constant_voltage")(time_step, constant_voltage);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_constant_power  (double const time_step, double const constant_power  )
{
    this->get_override("evolve_one_time_step_constant_power"  )(time_step, constant_power  );
}
void EnergyStorageDeviceWrap::evolve_one_time_step_constant_load   (double const time_step, double const constant_load   )
{
    this->get_override("evolve_one_time_step_constant_load"   )(time_step, constant_load   );
}
void EnergyStorageDeviceWrap::evolve_one_time_step_changing_current(double const time_step, double const changing_current)
{
    this->get_override("evolve_one_time_step_changing_current")(time_step, changing_current);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_changing_voltage(double const time_step, double const changing_voltage)
{
    this->get_override("evolve_one_time_step_changing_voltage")(time_step, changing_voltage);
}
void EnergyStorageDeviceWrap::evolve_one_time_step_changing_power  (double const time_step, double const changing_power  )
{
    this->get_override("evolve_one_time_step_changing_power"  )(time_step, changing_power  );
}
void EnergyStorageDeviceWrap::evolve_one_time_step_changing_load   (double const time_step, double const changing_load   )
{
    this->get_override("evolve_one_time_step_changing_load"   )(time_step, changing_load   );
}

double get_current(cap::EnergyStorageDevice const & dev)
{
    double current;
    dev.get_current(current);
    return current;
}

double get_voltage(cap::EnergyStorageDevice const & dev)
{
    double voltage;
    dev.get_voltage(voltage);
    return voltage;
}

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & python_object)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(python_object);
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(ptree);
    return cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));
}

// DEPRECATED
void ElectrochemicalImpedanceSpectroscopyData::impedance_spectroscopy(boost::python::object & python_device, boost::python::object & python_database)
{
    std::shared_ptr<cap::EnergyStorageDevice> device =
        boost::python::extract<std::shared_ptr<cap::EnergyStorageDevice>>(python_device);
    boost::property_tree::ptree const & database =
        boost::python::extract<boost::property_tree::ptree const &>(python_database);
    std::shared_ptr<boost::property_tree::ptree> eis_database =
        std::make_shared<boost::property_tree::ptree>(database.get_child("impedance_spectroscopy"));
#ifdef WITH_GSL
    std::map<double,std::complex<double>> eis_data = cap::impedance_spectroscopy(device, eis_database);
#else
    std::map<double,std::complex<double>> eis_data;
#endif
    data.insert(eis_data.begin(), eis_data.end());
}

void ElectrochemicalImpedanceSpectroscopyData::measure_impedance(boost::python::object & python_device, boost::python::object & python_database)
{
    std::shared_ptr<cap::EnergyStorageDevice> device =
        boost::python::extract<std::shared_ptr<cap::EnergyStorageDevice>>(python_device);
    boost::property_tree::ptree const & database =
        boost::python::extract<boost::property_tree::ptree const &>(python_database);
    std::shared_ptr<boost::property_tree::ptree> eis_database =
        std::make_shared<boost::property_tree::ptree>(database.get_child("impedance_spectroscopy"));
#ifdef WITH_GSL
    std::map<double,std::complex<double>> eis_data = cap::measure_impedance(device, eis_database);
#else
    std::map<double,std::complex<double>> eis_data;
#endif
    data.insert(eis_data.begin(), eis_data.end());
}

boost::python::list ElectrochemicalImpedanceSpectroscopyData::get_frequency() const
{
    boost::python::list frequency;
    for (auto const & p : data)
        frequency.append(p.first);
    return frequency;
}

boost::python::list ElectrochemicalImpedanceSpectroscopyData::get_complex_impedance() const
{
    boost::python::list complex_impedance;
    for (auto const & p : data)
        complex_impedance.append(p.second);
    return complex_impedance;
}

void ElectrochemicalImpedanceSpectroscopyData::clear()
{
    data.clear();
}

std::shared_ptr<boost::property_tree::ptree> compute_equivalent_circuit(boost::python::object & python_object)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(python_object);
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(ptree);
    std::shared_ptr<boost::property_tree::ptree> equivalent_circuit_database =
        std::make_shared<boost::property_tree::ptree>();
    cap::compute_equivalent_circuit(device_database, equivalent_circuit_database);
    return equivalent_circuit_database;
}
} // end namespace pycap

