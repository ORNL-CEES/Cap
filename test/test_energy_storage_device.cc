#define BOOST_TEST_MODULE EnergyStorageDevice
#define BOOST_TEST_MAIN

#include <cap/energy_storage_device.h>
#include <cap/resistor_capacitor.h>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>
#include <sstream>

//BOOST_SERIALIZATION_ASSUME_ABSTRACT(cap::EnergyStorageDevice)
//BOOST_CLASS_EXPORT(cap::SeriesRC)

// list of valid inputs to build an EnergyStorageDevice
// These are meant as example
std::list<std::string> const valid_device_input = {
    "series_rc.info",
    "parallel_rc.info",
#ifdef WITH_DEAL_II
    "super_capacitor.info",
#endif
    };

BOOST_AUTO_TEST_CASE( test_factory )
{
    for (auto const & filename : valid_device_input)
    {
        auto ptree = std::make_shared<boost::property_tree::ptree>();
        boost::property_tree::info_parser::read_info(filename, *ptree);
        BOOST_CHECK_NO_THROW(
            cap::buildEnergyStorageDevice(
                std::make_shared<cap::Parameters>(ptree) ) );
    }

    // invalid type must throw an exception
    auto ptree = std::make_shared<boost::property_tree::ptree>();
    ptree->put("type", "InvalidDeviceType");
    BOOST_CHECK_THROW(
        cap::buildEnergyStorageDevice(
            std::make_shared<cap::Parameters>(ptree) ),
        std::runtime_error );
}

// TODO: won't work for SuperCapacitor
#ifdef WITH_DEAL_II
    BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( test_serialization, 1 )
#endif

BOOST_AUTO_TEST_CASE( test_serialization )
{
    for (auto const & filename : valid_device_input)
    {
        auto ptree = std::make_shared<boost::property_tree::ptree>();
        boost::property_tree::info_parser::read_info(filename, *ptree);
        auto original_device = cap::buildEnergyStorageDevice(
            std::make_shared<cap::Parameters>(ptree) );

        original_device->evolve_one_time_step_constant_voltage(0.1, 2.1);
        double original_voltage;
        double original_current;
        original_device->get_voltage(original_voltage);
        original_device->get_current(original_current);

try {
        std::stringstream ss;
        // save device
        boost::archive::text_oarchive oa(ss);
        oa.register_type<cap::SeriesRC>();
        oa.register_type<cap::ParallelRC>();
        oa<<original_device;
        // print the content of the stream to the screen
        std::cout<<ss.str()<<"\n";
        BOOST_CHECK( !ss.str().empty() );

        // restore device
        boost::archive::text_iarchive ia(ss);
        ia.register_type<cap::SeriesRC>();
        ia.register_type<cap::ParallelRC>();
        std::shared_ptr<cap::EnergyStorageDevice> restored_device;
        ia>>restored_device;
        double restored_voltage;
        double restored_current;
        restored_device->get_voltage(restored_voltage);
        restored_device->get_current(restored_current);
        BOOST_CHECK_EQUAL(original_voltage, restored_voltage);
        BOOST_CHECK_EQUAL(original_current, restored_current);
} catch (boost::archive::archive_exception e) {
    BOOST_TEST_MESSAGE("unable to serialize the device");
    BOOST_TEST(false);
} catch (...) {
    throw std::runtime_error("unexpected exception occured");
}
    }

}
