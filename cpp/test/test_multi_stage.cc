#define BOOST_TEST_MODULE TestMultistage
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <memory>
#include <limits>
#include <cmath>

namespace cap {

class EndCriterion {
public:
    EndCriterion() { }
    virtual ~EndCriterion() { }
    virtual void reset(double const , std::shared_ptr<cap::EnergyStorageDevice> ) { };
    virtual bool check(double const time, std::shared_ptr<cap::EnergyStorageDevice> dev) const = 0;
};

class TimeLimit : public EndCriterion {
public:
    TimeLimit(std::shared_ptr<boost::property_tree::ptree const> database)
    {
        duration = database->get<double>("duration");
    }
    void reset(double const time, std::shared_ptr<cap::EnergyStorageDevice> )       override { tick = time; }
    bool check(double const time, std::shared_ptr<cap::EnergyStorageDevice> ) const override { return time - tick >= duration; }
    
private:
    double tick;
    double duration;
};

class VoltageLimit : public EndCriterion {
public:
    VoltageLimit(std::shared_ptr<boost::property_tree::ptree const> database, std::function<bool(double const, double const)> const & comp)
    : compare(comp)
    {
        voltage_limit = database->get<double>("voltage_limit");
    }
    bool check(double const , std::shared_ptr<cap::EnergyStorageDevice> dev) const override
    {
        double voltage;
        dev->get_voltage(voltage);
        return compare(voltage, voltage_limit);
    }
    double voltage_limit;
    std::function<bool(double const, double const)> compare;
};

class OperatingConditions {
public:
    OperatingConditions() { }
    virtual ~OperatingConditions() { }
    virtual void evolve_one_time_step(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev) = 0;
};

class Potentiostatic : public OperatingConditions {
public:
    Potentiostatic(std::shared_ptr<boost::property_tree::ptree const> database)
    {
        voltage = database->get<double>("voltage");
    }
    void evolve_one_time_step(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev) override
    {
        time += time_step;
        dev->evolve_one_time_step_constant_voltage(time_step, voltage);

    }
private:
    double voltage;
};

class Galvanostatic : public OperatingConditions {
public:
    Galvanostatic(std::shared_ptr<boost::property_tree::ptree const> database)
    {
        current = database->get<double>("current");
    }
    void evolve_one_time_step(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev) override
    {
        time += time_step;
        dev->evolve_one_time_step_constant_current(time_step, current);

    }
private:
    double current;
};

class Hold : public OperatingConditions {
public:
    Hold(std::shared_ptr<boost::property_tree::ptree const> )
    { }
    void evolve_one_time_step(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev) override
    {
        time += time_step;
        double voltage;
        dev->get_voltage(voltage);
        dev->evolve_one_time_step_constant_voltage(time_step, voltage);
    }
};

class Rest : public OperatingConditions {
public:
    Rest(std::shared_ptr<boost::property_tree::ptree const> )
    { }
    void evolve_one_time_step(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev) override
    {
        time += time_step;
        dev->evolve_one_time_step_constant_current(time_step, 0.0);
    }
};

class Stage {
public:
    Stage(std::shared_ptr<boost::property_tree::ptree const> database)
    {
        std::string type = database->get<std::string>("type");
        if (type.compare("potentiostatic") == 0)
            operating_conditions = std::make_shared<Potentiostatic>(database);
        else if (type.compare("galvanostatic" ) == 0)
            operating_conditions = std::make_shared<Galvanostatic >(database);
        else if (type.compare("hold"          ) == 0)
            operating_conditions = std::make_shared<Hold          >(database);
        else if (type.compare("rest"          ) == 0)
            operating_conditions = std::make_shared<Rest          >(database);
        else
            throw std::runtime_error("invalid operating conditions type "+type);

        type = database->get<std::string>("end_criterion");
        if (type.compare("time") == 0)
            end_criterion = std::make_shared<TimeLimit   >(database);
        else if (type.compare("voltage_greater_than") == 0)
            end_criterion = std::make_shared<VoltageLimit>(database, std::greater_equal<double>());
        else if (type.compare("voltage_less_than"   ) == 0)
            end_criterion = std::make_shared<VoltageLimit>(database, std::less_equal   <double>());
        else
            throw std::runtime_error("invalid end criterion "+type);
            
        max_steps    = database->get<int   >("max_steps"   , std::numeric_limits<int>::max());
        max_duration = database->get<double>("max_duration", std::nan("")                   );
    }
    virtual ~Stage() { }
    virtual int run(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev, std::ostream & os);
private:
    std::shared_ptr<OperatingConditions> operating_conditions;
    std::shared_ptr<EndCriterion       > end_criterion       ;
    double max_steps;
    double max_duration;
};

void report_data(double const time, std::shared_ptr<cap::EnergyStorageDevice const> dev, std::ostream & os)
{
    double current;
    double voltage;
    dev->get_current(current);
    dev->get_voltage(voltage);
    os<<boost::format("%20.15e  %20.15e  %20.15e  \n")
        % time
        % current
        % voltage
        ;
}

int
Stage::
run(double & time, double & time_step, std::shared_ptr<cap::EnergyStorageDevice> dev, std::ostream & os)
{
    double const tick = time;
    int    step = 0;
    end_criterion->reset(time, dev);
    while (true)
    {
        ++step;
        operating_conditions->evolve_one_time_step(time, time_step, dev);
        report_data(time, dev, os);
        if (end_criterion->check(time, dev))
            break;
        if (step >= max_steps)
            throw std::runtime_error("number of steps reach user specified maximum "+std::to_string(max_steps));
        if (time - tick >= max_duration) 
            throw std::runtime_error("duration exceed user specified maximum "+std::to_string(max_duration));
    }
    return step;
}

void operate(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os)
{
    int const stages = database->get<int>("stages");
    std::list<std::shared_ptr<cap::Stage> > list;
    for (int k = 0; k < stages; ++k) 
    {
        std::shared_ptr<boost::property_tree::ptree> stage_database =
            std::make_shared<boost::property_tree::ptree>(database->get_child("stage_"+std::to_string(k)));
        list.push_back(std::make_shared<cap::Stage>(stage_database));
    }


    double const initial_voltage = database->get<double>("initial_voltage");
    double       time_step       = database->get<double>("time_step"      );
    dev->reset_voltage(initial_voltage);
    double const initial_time = 0.0;
    double       time         = initial_time;
    int const cycles = database->get<int>("cycles");
    for (int n = 0; n < cycles; ++n)
    {
        double const cycle_tick = time;
        BOOST_FOREACH(std::shared_ptr<cap::Stage> stage, list)
            stage->run(time, time_step, dev, os);
        double const cycle_duration = time - cycle_tick;
        std::cout<<"cycle_"<<n<<" duration = "<<cycle_duration<<"\n";
    }
    double const final_time = time;
    double const duration = final_time - initial_time;
    std::cout<<"total duration = "<<duration<<"\n"; 
}


} // end namespace cap

BOOST_AUTO_TEST_CASE( test_multistage )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_multi_stage", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(boost::mpi::communicator(), *device_database);

    // operate the system
    std::fstream fout;
    fout.open("multi_stage_data", std::fstream::out);
    std::shared_ptr<boost::property_tree::ptree> operating_conditions_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("operating_conditions"));
    cap::operate(device, operating_conditions_database, fout);
    fout.close();
}
