#define BOOST_TEST_MODULE TestEquivalentCircuit
#define BOOST_TEST_MAIN
#include <cap/super_capacitor.h>
#include <cap/mp_values.h>
#include <cap/utils.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <functional>
#include <algorithm>


void report_time_current_voltage(std::vector<double> const & time,
    std::vector<double> const & current,
    std::vector<double> const & voltage,
    std::ostream & os = std::cout)
{
    std::size_t const n = time.size();
    if ((current.size() != n)
        || (voltage.size() != n))
        throw std::runtime_error("vector size mismatch");
    for (std::size_t i = 0; i < n; ++i)
        os<<boost::format("  %10.1f  %10.7f  %10.5f  \n")
            % time[i]
            % current[i]
            % voltage[i]
            ;
}

BOOST_AUTO_TEST_CASE( test_report_time_current_voltage )
{
    std::vector<double> t(1);
    std::vector<double> i(2);
    std::vector<double> u(1);
    BOOST_CHECK_THROW( report_time_current_voltage(t, i, u), std::runtime_error );
}

// forward declaration (implementation is in "common.cc")
std::shared_ptr<boost::property_tree::ptree> initialize_database();
void put_default_parameters(boost::property_tree::ptree & params);

BOOST_AUTO_TEST_CASE( test_equivalent_circuit )
{
    // create a database
    std::shared_ptr<boost::property_tree::ptree> database = initialize_database();
    database->put("verbose", false);
    std::shared_ptr<boost::property_tree::ptree> in(new boost::property_tree::ptree);
    put_default_parameters(*in);
    in->put("test_case",       2  );
    in->put("time_step",       0.1);
    in->put("initial_time",    0.0);
    in->put("final_time",     15.0);
    in->put("max_cycles",      1  );

{   // run high-fidelity model
    cap::SuperCapacitorProblem<2> high_fidelity_model(database);
    std::shared_ptr<boost::property_tree::ptree> out(new boost::property_tree::ptree);
    high_fidelity_model.run(in, out);
    
    std::vector<double> time    = cap::to_vector<double>(out->get<std::string>("time")   );
    std::vector<double> current = cap::to_vector<double>(out->get<std::string>("current"));
    std::vector<double> voltage = cap::to_vector<double>(out->get<std::string>("voltage"));
    report_time_current_voltage(time, current, voltage);
}



    database = in; // TODO: encapsulate 0-D model
    double const sandwich_height = database->get<double>("geometry.sandwich_height");
    double const cross_sectional_area = sandwich_height * 1.0;
    double const electrode_width = database->get<double>("geometry.electrode_width");
    double const separator_width = database->get<double>("geometry.separator_width");
    double const collector_width = database->get<double>("geometry.collector_width");

    // getting the material parameters values
    std::shared_ptr<boost::property_tree::ptree> material_properties_database = 
        std::make_shared<boost::property_tree::ptree>(database->get_child("material_properties"));
    std::shared_ptr<cap::SuperCapacitorMPValues<2> > mp_values = std::shared_ptr<cap::SuperCapacitorMPValues<2> >
        (new cap::SuperCapacitorMPValues<2>(cap::SuperCapacitorMPValuesParameters<2>(material_properties_database)));
    // build dummy cell itertor and set its material id
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube (triangulation);
    dealii::DoFHandler<2> dof_handler(triangulation);
    dealii::DoFHandler<2>::active_cell_iterator cell = 
        dof_handler.begin_active();
    // electrode
    cell->set_material_id(
        database->get<dealii::types::material_id>("material_properties.anode_electrode_material_id"));
    std::vector<double> electrode_solid_electrical_conductivity_values(1);
    std::vector<double> electrode_liquid_electrical_conductivity_values(1);
    mp_values->get_values("solid_electrical_conductivity", cell, electrode_solid_electrical_conductivity_values);
    mp_values->get_values("liquid_electrical_conductivity", cell, electrode_liquid_electrical_conductivity_values);
    double const electrode_resistivity = 
        1.0 / electrode_solid_electrical_conductivity_values[0]
        +
        1.0 / electrode_liquid_electrical_conductivity_values[0];
    double const electrode_resistance = electrode_resistivity * electrode_width / cross_sectional_area;
    std::vector<double> electrode_specific_capacitance_values(1);
    mp_values->get_values("specific_capacitance", cell, electrode_specific_capacitance_values);
    double const electrode_capacitance =
        electrode_specific_capacitance_values[0] * electrode_width * cross_sectional_area;
    // separator
    cell->set_material_id(
        database->get<dealii::types::material_id>("material_properties.separator_material_id"));
    std::vector<double> separator_liquid_electrical_conductivity_values(1);
    mp_values->get_values("liquid_electrical_conductivity", cell, separator_liquid_electrical_conductivity_values);
    double const separator_resistivity = 
        1.0 / separator_liquid_electrical_conductivity_values[0];
    double const separator_resistance = separator_resistivity * separator_width / cross_sectional_area;
    // collector
    cell->set_material_id(
        database->get<dealii::types::material_id>("material_properties.anode_collector_material_id"));
    std::vector<double> collector_solid_electrical_conductivity_values(1);
    mp_values->get_values("solid_electrical_conductivity", cell, collector_solid_electrical_conductivity_values);
    double const collector_resistivity = 
        1.0 / collector_solid_electrical_conductivity_values[0];
    double const collector_resistance = collector_resistivity * collector_width / cross_sectional_area;

    std::cout<<"electrode_capacitance="<<electrode_capacitance<<"\n";
    std::cout<<"electrode_resistance="<<electrode_resistance<<"\n";
    std::cout<<"separator_resistance="<<separator_resistance<<"\n";
    std::cout<<"collector_resistance="<<collector_resistance<<"\n";

    // compute the effective resistance and capacitance
//    these are the values I infered from the ``high-fidelity'' model
//    if I plug them into the 0-D model it produces the correct behavior
//    double const sandwich_capacitance = 0.02622;
//    double const sandwich_resistance = 59.96;
    double const sandwich_capacitance = electrode_capacitance / 2.0;
    double const sandwich_resistance = 2.0 * electrode_resistance + separator_resistance + 2.0 * collector_resistance;
    double const charge_current_density = database->get<double>("boundary_values.charge_current_density");
    double const discharge_current_density = database->get<double>("boundary_values.discharge_current_density");
    double const charge_current = charge_current_density * (collector_width * 1.0); // TODO: tab cross section
    double const discharge_current = discharge_current_density * (collector_width * 1.0); // TODO: 
    double const initial_voltage =  database->get<double>("boundary_values.initial_potential");
    std::cout<<"sandwich_capacitance="<<sandwich_capacitance<<"\n";
    std::cout<<"sandwich_resistance="<<sandwich_resistance<<"\n";
    std::cout<<"charge_current="<<charge_current<<"\n";
    std::cout<<"discharge_current="<<discharge_current<<"\n";
    std::cout<<"initial_voltage="<<initial_voltage<<"\n";
    

    // compute voltage evolution at constant current charge
    double const initial_time = database->get<double>("initial_time");
    double const final_time = database->get<double>("final_time");
    double const time_step = database->get<double>("time_step");
    std::size_t const n = std::ceil((final_time - initial_time) / time_step) + 1;
    std::vector<double> current(n, charge_current);
    std::vector<double> time(n);
    for (std::size_t i = 0; i < n; ++i)
        time[i] = initial_time + static_cast<double>(i) * (final_time - initial_time) / static_cast<double>(n-1) ;
    std::vector<double> voltage(n);
    cap::approximate_integral_with_trapezoidal_rule(time.begin(), time.end(), current.begin(), voltage.begin(), initial_voltage * sandwich_capacitance);
    std::transform(voltage.begin(), voltage.end(), current.begin(), voltage.begin(),
        [&sandwich_capacitance, &sandwich_resistance] (double const & U, double const & I)
            { return U / sandwich_capacitance + sandwich_resistance * I; }
        );
    report_time_current_voltage(time, current, voltage);

}

#include "common.cc"
