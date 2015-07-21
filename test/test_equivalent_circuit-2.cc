#define BOOST_TEST_MODULE TestEquivalentCircuit
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/mp_values.h>
#include <cap/utils.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <memory>
#include <cmath>

namespace cap {

// reads database for finite element model and write database for equivalent
// circuit model
void compute_equivalent_circuit(std::shared_ptr<boost::property_tree::ptree const> input_database,
                                std::shared_ptr<boost::property_tree::ptree      > output_database)
{
    // TODO: of course we could clear the database or just overwrite but for
    // now let's just throw an exception if it is not empty
    if (!output_database->empty())
        throw std::runtime_error("output_database was not empty...");

    double const sandwich_height = input_database->get<double>("geometry.sandwich_height");
    double const cross_sectional_area = sandwich_height * 1.0;
    double const electrode_width = input_database->get<double>("geometry.anode_electrode_width");
    double const separator_width = input_database->get<double>("geometry.separator_width");
    double const collector_width = input_database->get<double>("geometry.anode_collector_width");

    // getting the material parameters values
    std::shared_ptr<boost::property_tree::ptree> material_properties_database = 
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("material_properties"));
    cap::MPValuesParameters<2> mp_values_params(material_properties_database);
    std::shared_ptr<boost::property_tree::ptree> geometry_database = 
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("geometry"));
    mp_values_params.geometry = std::make_shared<cap::DummyGeometry<2>>(geometry_database);
    std::shared_ptr<cap::MPValues<2> > mp_values = std::shared_ptr<cap::MPValues<2> >
        (new cap::MPValues<2>(mp_values_params));
    // build dummy cell itertor and set its material id
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube (triangulation);
    dealii::DoFHandler<2> dof_handler(triangulation);
    dealii::DoFHandler<2>::active_cell_iterator cell = 
        dof_handler.begin_active();
    // electrode
    cell->set_material_id(
        input_database->get<dealii::types::material_id>("material_properties.anode_electrode_material_id"));
    std::vector<double> electrode_solid_electrical_conductivity_values(1);
    std::vector<double> electrode_liquid_electrical_conductivity_values(1);
    mp_values->get_values("solid_electrical_conductivity", cell, electrode_solid_electrical_conductivity_values);
    mp_values->get_values("liquid_electrical_conductivity", cell, electrode_liquid_electrical_conductivity_values);
    double const electrode_resistivity =
        ( 1.0 / electrode_solid_electrical_conductivity_values[0]
          +
          1.0 / electrode_liquid_electrical_conductivity_values[0] 
          +
          1.0 / ( electrode_solid_electrical_conductivity_values[0]
                  +
                  electrode_liquid_electrical_conductivity_values[0] )
        ) / 3.0
        ;
    double const electrode_resistance = electrode_resistivity * electrode_width / cross_sectional_area;
    std::vector<double> electrode_specific_capacitance_values(1);
    mp_values->get_values("specific_capacitance", cell, electrode_specific_capacitance_values);
    double const electrode_capacitance =
        electrode_specific_capacitance_values[0] * electrode_width * cross_sectional_area;
    std::vector<double> electrode_exchange_current_density_values(1);
    mp_values->get_values("faradaic_reaction_coefficient", cell, electrode_exchange_current_density_values);
    double const electrode_leakage_resistance = 1.0 / (electrode_exchange_current_density_values[0] * electrode_width * cross_sectional_area);
    std::cout<<"ELECTRODE\n";
    std::cout<<"    specific_capacitance="<<electrode_specific_capacitance_values[0]<<"\n";
    std::cout<<"    solid_electrical_conductivity="<<electrode_solid_electrical_conductivity_values[0]<<"\n";
    std::cout<<"    liquid_electrical_conductivity="<<electrode_liquid_electrical_conductivity_values[0]<<"\n";
    std::cout<<"    exchange_current_density="<<electrode_exchange_current_density_values[0]<<"\n";
    std::cout<<"    width="<<electrode_width<<"\n";
    std::cout<<"    cross_sectional_area="<<cross_sectional_area<<"\n";
    // separator
    cell->set_material_id(
        input_database->get<dealii::types::material_id>("material_properties.separator_material_id"));
    std::vector<double> separator_liquid_electrical_conductivity_values(1);
    mp_values->get_values("liquid_electrical_conductivity", cell, separator_liquid_electrical_conductivity_values);
    double const separator_resistivity = 
        1.0 / separator_liquid_electrical_conductivity_values[0];
    double const separator_resistance = separator_resistivity * separator_width / cross_sectional_area;
    std::cout<<"SEPARATOR\n";
    std::cout<<"    liquid_electrical_conductivity="<<separator_liquid_electrical_conductivity_values[0]<<"\n";
    std::cout<<"    width="<<separator_width<<"\n";
    std::cout<<"    cross_sectional_area="<<cross_sectional_area<<"\n";
    // collector
    cell->set_material_id(
        input_database->get<dealii::types::material_id>("material_properties.anode_collector_material_id"));
    std::vector<double> collector_solid_electrical_conductivity_values(1);
    mp_values->get_values("solid_electrical_conductivity", cell, collector_solid_electrical_conductivity_values);
    double const collector_resistivity = 
        1.0 / collector_solid_electrical_conductivity_values[0];
    double const collector_resistance = collector_resistivity * collector_width / cross_sectional_area;
    std::cout<<"COLLECTOR\n";
    std::cout<<"    solid_electrical_conductivity="<<collector_solid_electrical_conductivity_values[0]<<"\n";
    std::cout<<"    width="<<collector_width<<"\n";
    std::cout<<"    cross_sectional_area="<<cross_sectional_area<<"\n";

    std::cout<<"electrode_capacitance="<<electrode_capacitance<<"\n";
    std::cout<<"electrode_resistance="<<electrode_resistance<<"\n";
    std::cout<<"electrode_leakage_resistance="<<electrode_leakage_resistance<<"\n";
    std::cout<<"separator_resistance="<<separator_resistance<<"\n";
    std::cout<<"collector_resistance="<<collector_resistance<<"\n";

    // compute the effective resistance and capacitance
    double const sandwich_capacitance = electrode_capacitance / 2.0;
    double const sandwich_resistance = 2.0 * electrode_resistance + separator_resistance + 2.0 * collector_resistance;
    double const sandwich_leakage_resistance = 2.0 * electrode_leakage_resistance;
    std::cout<<"sandwich_capacitance="<<sandwich_capacitance<<"\n";
    std::cout<<"sandwich_resistance="<<sandwich_resistance<<"\n";
    std::cout<<"sandwich_leakage_resistance="<<sandwich_leakage_resistance<<"\n";

    output_database->put("capacitance"        , sandwich_capacitance       );
    output_database->put("series_resistance"  , sandwich_resistance        );
    output_database->put("parallel_resistance", sandwich_leakage_resistance);
    if (std::isfinite(sandwich_leakage_resistance))
        output_database->put("type", "ParallelRC");
    else
        output_database->put("type", "SeriesRC"  );
        
}



void report(double const time, std::shared_ptr<EnergyStorageDevice const> dev, std::ostream & os)
{
    double current;
    dev->get_current(current);
    double voltage;
    dev->get_voltage(voltage);
    os<<boost::format("  %22.15e  %22.15e  %22.15e  \n")
        % time
        % current
        % voltage
        ;
    // debug output of the solution fields and fluxes
    std::ostream null_sink(0);
    dev->print_data(null_sink);
}



void cycling_charge_discharge(std::shared_ptr<boost::property_tree::ptree const> database, std::shared_ptr<cap::EnergyStorageDevice> dev, std::ostream & os)
{
    double const time_step           = database->get<double>("time_step"               );
    double const initial_voltage     = database->get<double>("initial_voltage"         );
    double const charge_current      = database->get<double>("charge_current"          );
    double const discharge_current   = database->get<double>("discharge_current"       );
    double const voltage_upper_limit = database->get<double>("voltage_upper_limit"     );
    double const voltage_lower_limit = database->get<double>("voltage_lower_limit"     );
    double const charge_hold_time    = database->get<double>("charge_hold_time"   , 0.0);
    double const discharge_hold_time = database->get<double>("discharge_hold_time", 0.0);
    double const charge_rest_time    = database->get<double>("charge_rest_time"   , 0.0);
    double const discharge_rest_time = database->get<double>("discharge_rest_time", 0.0);
    int    const cycles              = database->get<int   >("cycles"                  );
    dev->reset_voltage(initial_voltage);
    double voltage = initial_voltage;
    double time = 0.0;
    double tick;
    for (int n = 0; n < cycles; ++n)
    {
        // glavanostatic charge
        while (voltage < voltage_upper_limit)
        {
            dev->evolve_one_time_step_constant_current(time_step, charge_current);
            dev->get_voltage(voltage);
            time += time_step;
            report(time, dev, os);
        }
        // potentiostatic hold
        tick = time;
        while (time - tick < charge_hold_time)
        {
            dev->evolve_one_time_step_constant_voltage(time_step, voltage);
            time += time_step;
            report(time, dev, os);
        }
        // rest at open circuit
        tick = time;
        while (time - tick < charge_rest_time)
        {
            dev->evolve_one_time_step_constant_current(time_step, 0.0);
            time += time_step;
            report(time, dev, os);
        }
        // glavanostatic discharge
        while (voltage > voltage_lower_limit)
        {
            dev->evolve_one_time_step_constant_current(time_step, discharge_current);
            dev->get_voltage(voltage);
            time += time_step;
            report(time, dev, os);
        }
        // potentiostatic hold
        tick = time;
        while (time - tick < discharge_hold_time)
        {
            dev->evolve_one_time_step_constant_voltage(time_step, voltage);
            time += time_step;
            report(time, dev, os);
        }
        // rest at open circuit
        tick = time;
        while (time - tick < discharge_rest_time)
        {
            dev->evolve_one_time_step_constant_current(time_step, 0.0);
            time += time_step;
            report(time, dev, os);
        }

    }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_equivalent_circuit )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    boost::property_tree::xml_parser::read_xml("input_equivalent_circuit", *input_database,
        boost::property_tree::xml_parser::trim_whitespace | boost::property_tree::xml_parser::no_comments);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // 
    std::shared_ptr<boost::property_tree::ptree> not_empty_database =
        std::make_shared<boost::property_tree::ptree>();
    not_empty_database->put("something", "not_empty");
    BOOST_CHECK_THROW(cap::compute_equivalent_circuit(device_database, not_empty_database), std::runtime_error);

    // build an equivalent resistor capacitor circuit
    std::shared_ptr<boost::property_tree::ptree> equivalent_circuit_database =
        std::make_shared<boost::property_tree::ptree>();
    cap::compute_equivalent_circuit(device_database, equivalent_circuit_database);
    std::shared_ptr<cap::EnergyStorageDevice> equivalent_circuit =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(equivalent_circuit_database));
    
    std::shared_ptr<boost::property_tree::ptree> cycling_charge_discharge_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("cycling_charge_discharge"));

    std::fstream fout;
    fout.open("equivalent_circuit_data_lf", std::fstream::out);
    cap::cycling_charge_discharge(cycling_charge_discharge_database, equivalent_circuit, fout);

    fout.close();
    fout.open("equivalent_circuit_data_hf", std::fstream::out);
    cap::cycling_charge_discharge(cycling_charge_discharge_database, device            , fout);
    
}
