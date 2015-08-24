#include <cap/energy_storage_device.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/python.hpp>
#include <string>
#include <memory>


namespace cap {

struct EnergyStorageDeviceWrap: EnergyStorageDevice, boost::python::wrapper<EnergyStorageDevice>
{
    void get_voltage(double & voltage) const
    {
        this->get_override("get_voltage")(voltage);
    }
    void get_current(double & current) const
    {
        this->get_override("get_current")(current);
    }
    void evolve_one_time_step_constant_current(double const time_step, double const constant_current)
    {
        this->get_override("evolve_one_time_step_constant_current")(time_step, constant_current);
    }
    void evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage)
    {
        this->get_override("evolve_one_time_step_constant_voltage")(time_step, constant_voltage);
    }
    void evolve_one_time_step_constant_power  (double const time_step, double const constant_power  )
    {
        this->get_override("evolve_one_time_step_constant_power"  )(time_step, constant_power  );
    }
    void evolve_one_time_step_constant_load   (double const time_step, double const constant_load   )
    {
        this->get_override("evolve_one_time_step_constant_load"   )(time_step, constant_load   );
    }
    void evolve_one_time_step_changing_current(double const time_step, double const changing_current)
    {
        this->get_override("evolve_one_time_step_changing_current")(time_step, changing_current);
    }
    void evolve_one_time_step_changing_voltage(double const time_step, double const changing_voltage)
    {
        this->get_override("evolve_one_time_step_changing_voltage")(time_step, changing_voltage);
    }
    void evolve_one_time_step_changing_power  (double const time_step, double const changing_power  )
    {
        this->get_override("evolve_one_time_step_changing_power"  )(time_step, changing_power  );
    }
    void evolve_one_time_step_changing_load   (double const time_step, double const changing_load   )
    {
        this->get_override("evolve_one_time_step_changing_load"   )(time_step, changing_load   );
    }
};

double get_current(EnergyStorageDevice const & dev)
{
    double current;
    dev.get_current(current);
    return current;
}

double get_voltage(EnergyStorageDevice const & dev)
{
    double voltage;
    dev.get_voltage(voltage);
    return voltage;
}

std::shared_ptr<EnergyStorageDevice>
build_energy_storage_device(std::string const & filename)
{
    boost::property_tree::ptree input_database;
    boost::property_tree::xml_parser::read_xml(filename, input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database.get_child("device"));
    return buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));
}

} // end namespace cap

BOOST_PYTHON_MODULE(pycap)
{
    boost::python::class_<cap::EnergyStorageDeviceWrap, std::shared_ptr<cap::EnergyStorageDeviceWrap>, boost::noncopyable>("EnergyStorageDevice", boost::python::no_init)
        .def("__init__", boost::python::make_constructor(&cap::build_energy_storage_device) )
        .def("get_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::get_voltage) )
        .def("get_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::get_current) )
        .def("evolve_one_time_step_constant_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_current) )
        .def("evolve_one_time_step_constant_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_voltage) )
        .def("evolve_one_time_step_constant_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_power  ) )
        .def("evolve_one_time_step_constant_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_constant_load   ) )
        .def("evolve_one_time_step_changing_current", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_current) )
        .def("evolve_one_time_step_changing_voltage", boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_voltage) )
        .def("evolve_one_time_step_changing_power"  , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_power  ) )
        .def("evolve_one_time_step_changing_load"   , boost::python::pure_virtual(&cap::EnergyStorageDevice::evolve_one_time_step_changing_load   ) )
        ;
    boost::python::register_ptr_to_python<std::shared_ptr<cap::EnergyStorageDevice>>();

    boost::python::def("get_current", cap::get_current);
    boost::python::def("get_voltage", cap::get_voltage);
}
