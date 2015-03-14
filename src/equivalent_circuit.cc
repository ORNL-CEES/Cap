#include <cap/equivalent_circuit.h>
#include <cap/utils.h>
#include <string>

namespace cap {

EquivalentCircuit::
EquivalentCircuit(std::shared_ptr<boost::property_tree::ptree const> )
{ }


void
EquivalentCircuit::
process_solution(double * data)
{
    data[JOULE_HEATING] = equivalent_circuit->R * equivalent_circuit->I * equivalent_circuit->I;
    data[VOLTAGE] = equivalent_circuit->U;
    data[CURRENT] = equivalent_circuit->I;
}



void
EquivalentCircuit::
evolve_one_time_step(double const & delta_t)
{
    if (capacitor_state == cap::GalvanostaticCharge)
        equivalent_circuit->evolve_one_time_step_constant_current(delta_t, I_charge);
    else if (capacitor_state == cap::GalvanostaticDischarge)
        equivalent_circuit->evolve_one_time_step_constant_current(delta_t, I_discharge);
    else if (capacitor_state == cap::PotentiostaticCharge)
        equivalent_circuit->evolve_one_time_step_constant_voltage(delta_t, U_charge);
    else if (capacitor_state == cap::PotentiostaticDischarge)
        equivalent_circuit->evolve_one_time_step_constant_voltage(delta_t, U_discharge);
    else if (capacitor_state == cap::Relaxation)
        equivalent_circuit->evolve_one_time_step_constant_current(delta_t, 0.0);
    else
        throw std::runtime_error("invalid capacitor state");
}



void
EquivalentCircuit::
reset(std::shared_ptr<boost::property_tree::ptree const> database)
{
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
//    // parallel
//    double const electrode_resistivity = 
//        1.0 / 
//        ( electrode_solid_electrical_conductivity_values[0]
//          +
//          electrode_liquid_electrical_conductivity_values[0] )
//        ;
    // series
    double const electrode_resistivity =
        ( 1.0 / electrode_solid_electrical_conductivity_values[0]
          +
          1.0 / electrode_liquid_electrical_conductivity_values[0] )
        * 0.33 // TODO:
        ;
    double const electrode_resistance = electrode_resistivity * electrode_width / cross_sectional_area;
    std::vector<double> electrode_specific_capacitance_values(1);
    mp_values->get_values("specific_capacitance", cell, electrode_specific_capacitance_values);
    double const electrode_capacitance =
        electrode_specific_capacitance_values[0] * electrode_width * cross_sectional_area;
    std::cout<<"ELECTRODE\n";
    std::cout<<"    specific_capacitance="<<electrode_specific_capacitance_values[0]<<"\n";
    std::cout<<"    solid_electrical_conductivity="<<electrode_solid_electrical_conductivity_values[0]<<"\n";
    std::cout<<"    liquid_electrical_conductivity="<<electrode_liquid_electrical_conductivity_values[0]<<"\n";
    std::cout<<"    width="<<electrode_width<<"\n";
    std::cout<<"    cross_sectional_area="<<cross_sectional_area<<"\n";
    // separator
    cell->set_material_id(
        database->get<dealii::types::material_id>("material_properties.separator_material_id"));
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
        database->get<dealii::types::material_id>("material_properties.anode_collector_material_id"));
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

    std::shared_ptr<boost::property_tree::ptree> db =
        std::make_shared<boost::property_tree::ptree>();
    db->put("capacitance"      , sandwich_capacitance);
    db->put("series_resistance", sandwich_resistance );

    equivalent_circuit = 
        std::make_shared<cap::SeriesRC>(
            std::make_shared<cap::Parameters>(db)
        );
    equivalent_circuit->reset_voltage(initial_voltage);
    I_charge    = charge_current;
    I_discharge = discharge_current;
}



void
EquivalentCircuit::
setup(CapacitorState const change)
{
    capacitor_state = change;
    if (capacitor_state == GalvanostaticCharge)
        equivalent_circuit->I = I_charge;
    else if (capacitor_state == GalvanostaticDischarge)
        equivalent_circuit->I = I_discharge;
    else if (capacitor_state == Relaxation)
        equivalent_circuit->I = 0.0;
    else
        throw std::runtime_error("invalid capacitor state");
    equivalent_circuit->reset(equivalent_circuit->U_C);
}



void 
EquivalentCircuit::
run(std::shared_ptr<boost::property_tree::ptree const> input_params,
    std::shared_ptr<boost::property_tree::ptree> output_params)
{
    this->reset(input_params);
    double       time_step    = input_params->get<double>("time_step"   );
    double const initial_time = input_params->get<double>("initial_time");
    double const final_time   = input_params->get<double>("final_time"  );
    std::size_t const max_cycles = input_params->get<std::size_t>("max_cycles");

    double const max_voltage = input_params->get<double>("boundary_values.charge_potential");
    double const min_voltage = input_params->get<double>("boundary_values.discharge_potential");

    std::vector<double> heat_production;
    std::vector<double> voltage;
    std::vector<double> current;
    std::vector<double> time;
    std::vector<std::string> capacitor_state;
    std::vector<int> cycle;

    double current_time = initial_time;
    std::size_t step = 0;
    std::size_t current_cycle = 0;
    double data[N_DATA];

    this->process_solution(data);
    heat_production.push_back(data[JOULE_HEATING]);
    voltage.push_back(data[VOLTAGE]);
    current.push_back(data[CURRENT]);
    time.push_back(current_time);
    cycle.push_back(current_cycle);
    capacitor_state.push_back("initialize");

    while (current_cycle < max_cycles) {
        ++current_cycle;
    
        this->setup(GalvanostaticCharge);
        while (current_time < final_time)
        {
            ++step;
            current_time += time_step;
            this->evolve_one_time_step(time_step);
            this->process_solution(data);
            heat_production.push_back(data[JOULE_HEATING]);
            voltage.push_back(data[VOLTAGE]);
            current.push_back(data[CURRENT]);
            time.push_back(current_time);
            cycle.push_back(current_cycle);
            capacitor_state.push_back("charging");
            if (data[VOLTAGE] >= max_voltage)
                break;
        }

        this->setup(GalvanostaticDischarge);
        while (current_time < final_time)
        {
            ++step;
            current_time += time_step;
            this->evolve_one_time_step(time_step);
            this->process_solution(data);
            heat_production.push_back(data[JOULE_HEATING]);
            voltage.push_back(data[VOLTAGE]);
            current.push_back(data[CURRENT]);
            time.push_back(current_time);
            cycle.push_back(current_cycle);
            capacitor_state.push_back("discharging");
            if (data[VOLTAGE] <= min_voltage)
                break;
        }
    }

    output_params->put("heat_production", to_string(heat_production));
    output_params->put("voltage",         to_string(voltage)        );
    output_params->put("current",         to_string(current)        );
    output_params->put("time",            to_string(time)           );
    output_params->put("capacitor_state", to_string(capacitor_state));
    output_params->put("cycle",           to_string(cycle)          );

}

} // end namespace cap

