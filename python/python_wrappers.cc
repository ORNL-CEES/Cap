#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/python.hpp>
#include <string>
#include <memory>


namespace cap {

struct EnergyStorageDeviceWrap : EnergyStorageDevice, boost::python::wrapper<EnergyStorageDevice>
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

double get_double(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return ptree.get<double>(path);
}

std::string get_string(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return ptree.get<std::string>(path);
}

int get_int(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return ptree.get<int>(path);
}

bool get_bool(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return ptree.get<bool>(path);
}

void put_double(boost::property_tree::ptree & ptree, std::string const & path, double const & value)
{
    ptree.put(path, value);
}

void put_string(boost::property_tree::ptree & ptree, std::string const & path, std::string const & value)
{
    ptree.put(path, value);
}

void put_int(boost::property_tree::ptree & ptree, std::string const & path, int const & value)
{
    ptree.put(path, value);
}

void put_bool(boost::property_tree::ptree & ptree, std::string const & path, bool const & value)
{
    ptree.put(path, value);
}

void parse_xml(boost::property_tree::ptree & ptree, std::string const & filename)
{
    boost::property_tree::xml_parser::read_xml(filename, ptree,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);
}

boost::property_tree::ptree get_child(boost::property_tree::ptree & ptree, std::string const & path)
{
    return ptree.get_child(path);
}

std::shared_ptr<EnergyStorageDevice>
build_energy_storage_device(boost::python::object & python_object)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(python_object);
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(ptree.get_child("device"));
    return buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));
}

struct ElectrochemicalImpedanceSpectroscopyData {
    ElectrochemicalImpedanceSpectroscopyData(std::map<double,std::complex<double>> const & d)
        : data(d)
    { }
    boost::python::list get_frequencies()
    {
        boost::python::list frequencies;
        for (auto p : data)
            frequencies.append(p.first);
        return frequencies;
    }
    boost::python::list get_complex_impedance()
    {
        boost::python::list frequencies;
        for (auto p : data)
            frequencies.append(p.second);
        return frequencies;
    }
    std::map<double,std::complex<double>> data;
};

std::shared_ptr<ElectrochemicalImpedanceSpectroscopyData>
measure_complex_impedance(boost::python::object & python_device, boost::python::object & python_database)
{
    std::ignore = python_database;
    std::shared_ptr<cap::EnergyStorageDevice> device =
        boost::python::extract<std::shared_ptr<cap::EnergyStorageDevice>>(python_device);
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(python_database);
    std::shared_ptr<boost::property_tree::ptree> eis_database =
        std::make_shared<boost::property_tree::ptree>(ptree.get_child("impedance_spectroscopy"));
    return std::make_shared<ElectrochemicalImpedanceSpectroscopyData>(impedance_spectroscopy(device, eis_database));
}


} // end namespace cap

BOOST_PYTHON_MODULE(pycap)
{
    boost::python::class_<cap::ElectrochemicalImpedanceSpectroscopyData, std::shared_ptr<cap::ElectrochemicalImpedanceSpectroscopyData>, boost::noncopyable>("ElectrochemicalImpedanceSpectroscopyData", boost::python::no_init)
        .def("__init__", boost::python::make_constructor(&cap::measure_complex_impedance) )
        .def("get_frequencies", &cap::ElectrochemicalImpedanceSpectroscopyData::get_frequencies)
        .def("get_complex_impedance", &cap::ElectrochemicalImpedanceSpectroscopyData::get_complex_impedance)
        ;
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

    boost::python::class_<boost::property_tree::ptree, std::shared_ptr<boost::property_tree::ptree>>("PropertyTree")
        .def("get_double", &cap::get_double)
        .def("get_string", &cap::get_string)
        .def("get_int"   , &cap::get_int   )
        .def("get_bool"  , &cap::get_bool  )
        .def("put_double", &cap::put_double)
        .def("put_string", &cap::put_string)
        .def("put_int"   , &cap::put_int   )
        .def("put_bool"  , &cap::put_bool  )
        .def("parse_xml" , &cap::parse_xml )
        .def("get_child" , &cap::get_child )
        ;
}
