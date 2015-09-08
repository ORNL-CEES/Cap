#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <cap/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/python.hpp>
#include <string>
#include <memory>
#include <map>


namespace pycap {

struct EnergyStorageDeviceWrap : cap::EnergyStorageDevice, boost::python::wrapper<cap::EnergyStorageDevice>
{
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

template <typename T>
boost::python::list get_array(boost::property_tree::ptree const & ptree, std::string const & path)
{
    std::vector<T> vector = 
        cap::to_vector<T>(ptree.get<std::string>(path));
    boost::python::list list;
    for (auto const & x : vector)
        list.append(x);
    return list;
}

boost::python::list get_array_double(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return get_array<double>(ptree, path);
}

boost::python::list get_array_string(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return get_array<std::string>(ptree, path);
}

boost::python::list get_array_int(boost::property_tree::ptree const & ptree, std::string const & path)
{
    return get_array<int>(ptree, path);
}

//boost::python::list get_array_bool(boost::property_tree::ptree const & ptree, std::string const & path)
//{
//    return get_array<bool>(ptree, path);
//}
//
void parse_xml(boost::property_tree::ptree & ptree, std::string const & filename)
{
    boost::property_tree::xml_parser::read_xml(filename, ptree,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);
}

void parse_json(boost::property_tree::ptree & ptree, std::string const & filename)
{
    boost::property_tree::json_parser::read_json(filename, ptree);
}

boost::property_tree::ptree get_child(boost::property_tree::ptree & ptree, std::string const & path)
{
    return ptree.get_child(path);
}

std::shared_ptr<cap::EnergyStorageDevice>
build_energy_storage_device(boost::python::object & python_object)
{
    boost::property_tree::ptree const & ptree =
        boost::python::extract<boost::property_tree::ptree const &>(python_object);
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(ptree.get_child("device"));
    return cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));
}

struct ElectrochemicalImpedanceSpectroscopyData {
    void impedance_spectroscopy(boost::python::object & python_device, boost::python::object & python_database)
    {
        std::shared_ptr<cap::EnergyStorageDevice> device =
            boost::python::extract<std::shared_ptr<cap::EnergyStorageDevice>>(python_device);
        boost::property_tree::ptree const & database =
            boost::python::extract<boost::property_tree::ptree const &>(python_database);
        std::shared_ptr<boost::property_tree::ptree> eis_database =
            std::make_shared<boost::property_tree::ptree>(database.get_child("impedance_spectroscopy"));
        std::map<double,std::complex<double>> eis_data = cap::impedance_spectroscopy(device, eis_database);
        data.insert(eis_data.begin(), eis_data.end());
    }
    void measure_impedance(boost::python::object & python_device, boost::python::object & python_database)
    {
        std::shared_ptr<cap::EnergyStorageDevice> device =
            boost::python::extract<std::shared_ptr<cap::EnergyStorageDevice>>(python_device);
        boost::property_tree::ptree const & database =
            boost::python::extract<boost::property_tree::ptree const &>(python_database);
        std::shared_ptr<boost::property_tree::ptree> eis_database =
            std::make_shared<boost::property_tree::ptree>(database.get_child("impedance_spectroscopy"));
        std::map<double,std::complex<double>> eis_data = cap::measure_impedance(device, eis_database);
        data.insert(eis_data.begin(), eis_data.end());
    }
    boost::python::list get_frequency() const
    {
        boost::python::list frequency;
        for (auto const & p : data)
            frequency.append(p.first);
        return frequency;
    }
    boost::python::list get_complex_impedance() const
    {
        boost::python::list complex_impedance;
        for (auto const & p : data)
            complex_impedance.append(p.second);
        return complex_impedance;
    }
    void clear()
    {
        data.clear();
    }
    std::map<double,std::complex<double>> data;
};

} // end namespace pycap

BOOST_PYTHON_MODULE(pycap)
{
    boost::python::class_<pycap::ElectrochemicalImpedanceSpectroscopyData, std::shared_ptr<pycap::ElectrochemicalImpedanceSpectroscopyData>>("ElectrochemicalImpedanceSpectroscopyData")
        .def("impedance_spectroscopy", &pycap::ElectrochemicalImpedanceSpectroscopyData::impedance_spectroscopy)
        .def("measure_impedance", &pycap::ElectrochemicalImpedanceSpectroscopyData::measure_impedance)
        .def("get_frequency", &pycap::ElectrochemicalImpedanceSpectroscopyData::get_frequency)
        .def("get_complex_impedance", &pycap::ElectrochemicalImpedanceSpectroscopyData::get_complex_impedance)
        .def("clear", &pycap::ElectrochemicalImpedanceSpectroscopyData::clear)
        ;
    boost::python::class_<pycap::EnergyStorageDeviceWrap, std::shared_ptr<pycap::EnergyStorageDeviceWrap>, boost::noncopyable>("EnergyStorageDevice", boost::python::no_init)
        .def("__init__", boost::python::make_constructor(&pycap::build_energy_storage_device) )
        .def("get_voltage", (&pycap::get_voltage) )
        .def("get_current", (&pycap::get_current) )
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

    boost::python::def("get_current", pycap::get_current);
    boost::python::def("get_voltage", pycap::get_voltage);

    boost::python::class_<boost::property_tree::ptree, std::shared_ptr<boost::property_tree::ptree>>("PropertyTree")
        .def("get_double", &pycap::get_double)
        .def("get_string", &pycap::get_string)
        .def("get_int"   , &pycap::get_int   )
        .def("get_bool"  , &pycap::get_bool  )
        .def("put_double", &pycap::put_double)
        .def("put_string", &pycap::put_string)
        .def("put_int"   , &pycap::put_int   )
        .def("put_bool"  , &pycap::put_bool  )
        .def("parse_xml" , &pycap::parse_xml )
        .def("parse_json", &pycap::parse_json)
        .def("get_child" , &pycap::get_child )
        .def("get_array_double", &pycap::get_array_double)
        .def("get_array_string", &pycap::get_array_string)
        .def("get_array_int"   , &pycap::get_array_int   )
//        .def("get_array_bool"  , &pycap::get_array_bool  )
        ;
}
