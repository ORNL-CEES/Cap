#include <cap/super_capacitor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <boost/format.hpp>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <functional>

namespace cap {

template <int dim>                         
SuperCapacitorProblem<dim>::               
SuperCapacitorProblem(std::shared_ptr<boost::property_tree::ptree const> database)
    : dof_handler(triangulation)
    , verbose(database->get<bool>("verbose", false))
    , symmetric_correction(true)
{ 
//    std::cout<<"initialize...\n";
    this->build_triangulation(database);
    this->set_cell_material_ids_and_face_boundary_ids(database);
    this->initialize_system(database);
}

template <int dim>
void
SuperCapacitorProblem<dim>::
run(std::shared_ptr<boost::property_tree::ptree const> input_params,
    std::shared_ptr<boost::property_tree::ptree>       output_params)
{                                               
    int const test_case = input_params->get<int>("test_case");
    if (test_case == 1) {
        this->run_constant_current_charge_constant_voltage_discharge(input_params, output_params);
    } else if (test_case == 2) {
        this->run_constant_current_cycling                          (input_params, output_params);
    } else {
       std::runtime_error("Unrecognized test case");
    } // end if
}

template <int dim>
void
SuperCapacitorProblem<dim>::
run_constant_current_cycling
    ( std::shared_ptr<boost::property_tree::ptree const> input_params
    , std::shared_ptr<boost::property_tree::ptree>       output_params
    )
{                                               
//    std::cout<<"run...\n";                      
    this->reset(input_params);

    unsigned int const thermal_block         = 1;
    unsigned int const electrochemical_block = 0;
    dealii::Vector<double> & thermal_solution         = this->solution.block(thermal_block        );
    dealii::Vector<double> & electrochemical_solution = this->solution.block(electrochemical_block);

   

{   // find initial solution for electrochemical
    electrochemical_solution = 0.0;
    double const dummy_time_step = 1.0;
    this->electrochemical_setup_system(dummy_time_step, Initialize);
    unsigned int step = 0;
    double solution_norm;
    double old_solution_norm = 0.0;
    while (true) {
        ++step;
        this->electrochemical_evolve_one_time_step(dummy_time_step);
        solution_norm = electrochemical_solution.l2_norm();
        if (std::abs(solution_norm - old_solution_norm) < 1.0e-8) {
            break;
        } // end if
        old_solution_norm = solution_norm;
    } // end while
    if (this->verbose)
        std::cout<<step<<" iterations for finding initial electrochemical solution\n";
}

    thermal_solution = 0.0; // TODO: initialize to ambient temperature
    
    double const max_voltage = input_params->get<double>("boundary_values.charge_potential");
    double const min_voltage = input_params->get<double>("boundary_values.discharge_potential");

    double const time_step    = input_params->get<double>("time_step"   );
    double const initial_time = input_params->get<double>("initial_time");
    double const final_time   = input_params->get<double>("final_time"  );
    unsigned int const max_cycles = input_params->get<unsigned int>("max_cycles");

    std::vector<double> max_temperature;
    std::vector<double> heat_production;
    std::vector<double> voltage;
    std::vector<double> current;
    std::vector<double> time;
    std::vector<std::string> capacitor_state;
    std::vector<int> cycle;
    double volume;
    double mass;

    double current_time = initial_time;
    unsigned int step          = 0;
    unsigned int current_cycle = 0;
    double data[N_DATA];

    this->thermal_setup_system(time_step);

    this->process_solution(data);
    this->report_data(current_time, data);
    max_temperature.push_back(data[TEMPERATURE]);
    heat_production.push_back(data[JOULE_HEATING]);
    voltage.push_back(data[VOLTAGE]);
    current.push_back(data[CURRENT]);
    time.push_back(current_time);
    cycle.push_back(current_cycle);
    capacitor_state.push_back("initialize");

    while (current_cycle < max_cycles) {
        if (current_time > final_time)
            break;
        ++current_cycle;
    
        this->electrochemical_setup_system(time_step, GalvanostaticCharge);
        while (current_time <= final_time) {
            ++step;
            current_time += time_step;
            this->electrochemical_evolve_one_time_step(time_step);
            this->thermal_evolve_one_time_step        (time_step);
            this->process_solution(data);
            this->report_data(current_time, data);
            max_temperature.push_back(data[TEMPERATURE]);
            heat_production.push_back(data[JOULE_HEATING]);
            voltage.push_back(data[VOLTAGE]);
            current.push_back(data[CURRENT]);
            time.push_back(current_time);
            cycle.push_back(current_cycle);
            capacitor_state.push_back("charging");
            volume = data[VOLUME];
            mass   = data[MASS  ];
            if (data[VOLTAGE] >= max_voltage) {
                break;
            } // end if
        } // end while

        this->electrochemical_setup_system(time_step, GalvanostaticDischarge);
        while (current_time <= final_time) {
            ++step;
            current_time += time_step;
            this->electrochemical_evolve_one_time_step(time_step);
            this->thermal_evolve_one_time_step        (time_step);
            this->process_solution(data);
            this->report_data(current_time, data);
            max_temperature.push_back(data[TEMPERATURE]);
            heat_production.push_back(data[JOULE_HEATING]);
            voltage.push_back(data[VOLTAGE]);
            current.push_back(data[CURRENT]);
            time.push_back(current_time);
            cycle.push_back(current_cycle);
            capacitor_state.push_back("discharging");
            if (data[VOLTAGE] <= min_voltage) {
                break;
            } // end if
        } // end while

    } // end while

    output_params->put("max_temperature", to_string(max_temperature));
    output_params->put("heat_production", to_string(heat_production));
    output_params->put("voltage",         to_string(voltage)        );
    output_params->put("current",         to_string(current)        );
    output_params->put("time",            to_string(time)           );
    output_params->put("capacitor_state", to_string(capacitor_state));
    output_params->put("cycle",           to_string(cycle)          );
    output_params->put("volume",          volume                    );
    output_params->put("mass",            mass                      );

    std::size_t const n = time.size();
    // compute power
    std::vector<double> power(n);
    std::transform(voltage.begin(), voltage.end(), current.begin(), power.begin(), 
        [&mass](double const & U, double const & I) { return U * I / mass; });
    // compute energy
    std::vector<double> energy(n);
    energy[0] = 0.0;
    for (std::size_t i = 1; i < n; ++i)
        if (capacitor_state[i].compare(capacitor_state[i-1]) == 0)
            energy[i] = energy[i-1] + 0.5 * (time[i] - time[i-1]) * (power[i] + power[i-1]);
        else
            energy[i] = 0.0; // reset to zero when toggle from charge to discharge
    // compute thermal losses
    std::vector<double> thermal_energy_loss(n);
    thermal_energy_loss[0] = 0.0;
    for (std::size_t i = 1; i < n; ++i)
        thermal_energy_loss[i] = thermal_energy_loss[i-1] + 0.5 * (time[i] - time[i-1]) * (heat_production[i] + heat_production[i-1]);
    // compute efficiency
    std::vector<double> efficiency(n);
    std::transform(power.begin(), power.end(), heat_production.begin(), efficiency.begin(), 
        [](double const & P, double const & Q) { return 100.0 * (std::abs(P) - Q) / std::abs(P); });

    output_params->put("power",               to_string(power)              );
    output_params->put("energy",              to_string(energy)             );
    output_params->put("thermal_energy_loss", to_string(thermal_energy_loss));
    output_params->put("efficiency",          to_string(efficiency)         );
    
    auto minmax_energy = std::minmax_element(energy.begin(), energy.end());
    double const seconds_per_hour = 3600.0;
    double const qoi_energy_density = 
         - *minmax_energy.first / seconds_per_hour;

    double const qoi_power_density = 
        std::accumulate(power.begin(), power.end(), 0.0, 
            [](double const & sum, double const & val) { return sum + std::abs(val); }
        ) / static_cast<double>(power.size());

    double const qoi_max_temperature = *std::max_element(max_temperature.begin(), max_temperature.end());
    double const qoi_efficiency =
          - *minmax_energy.first / *minmax_energy.second;

    output_params->put("quantities_of_interest.max_temperature", qoi_max_temperature);
    output_params->put("quantities_of_interest.energy_density",  qoi_energy_density );
    output_params->put("quantities_of_interest.power_density",   qoi_power_density  );
    output_params->put("quantities_of_interest.efficiency",      qoi_efficiency     );

    double max_heat_production = *std::max_element(heat_production.begin(), heat_production.end());
    std::vector<double> dummy(n, 0.0);
    for (std::size_t i = 1; i < n; ++i)
        dummy[i] = dummy[i-1] + (time[i] - time[i-1]) * max_heat_production;
    std::ofstream fout;
    fout.open("output_test_constant_current_cycling");
    fout<<boost::format("# %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  \n")
        % "time[s]"
        % "capacitor_state"
        % "cycle"
        % "current[A]"
        % "voltage[V]"
        % "max_temperature[K]"
        % "heat_production[W]"
        % "power[W/kg]"
        % "efficiency[%]"
        % "energy[Wh/kg]"
        % "loss[Wh/kg]"
        % "dummy[Wh/kg]"
    ;
    for (std::size_t i = 0; i < n; ++i)
          fout<<boost::format("  %7.1f  %15s  %5d  %10f  %10.3f  %18.3f  %18.4e  %11.4e  %13.1f  %13.4e  %11.4e  %12.4e  \n")
              % time[i]
              % capacitor_state[i]
              % cycle[i]
              % current[i]
              % voltage[i]
              % max_temperature[i]
              % heat_production[i]
              % power[i]
              % efficiency[i]
              % energy[i]
              % thermal_energy_loss[i]
              % dummy[i]
          ;
    fout.close();
}                                               

template <int dim>
void
SuperCapacitorProblem<dim>::
report_data(double time, double const * data)
{
    if (this->verbose) {
        std::cout<<boost::format("t=%6.1fs  U=%6.4fV  I=%6.4fA  Q=%8.6eW  T=%6.3fK  \n")
            % time
            % data[VOLTAGE]
            % data[CURRENT]
            % data[JOULE_HEATING]
            % data[TEMPERATURE]
        ;
    }
}

template <int dim>
void
SuperCapacitorProblem<dim>::
run_constant_current_charge_constant_voltage_discharge
    ( std::shared_ptr<boost::property_tree::ptree const> input_params
    , std::shared_ptr<boost::property_tree::ptree>       output_params
    )
{                                               
    std::cout<<"run...\n";                      
    this->reset(input_params);

    unsigned int const thermal_block         = 1;
    unsigned int const electrochemical_block = 0;
    dealii::Vector<double> & thermal_solution         = this->solution.block(thermal_block        );
    dealii::Vector<double> & electrochemical_solution = this->solution.block(electrochemical_block);

{   // find initial solution for electrochemical
    electrochemical_solution = 0.0;
    double const dummy_time_step = 1.0;
    this->electrochemical_setup_system(dummy_time_step, Initialize);
    unsigned int step = 0;
    double solution_norm;
    double old_solution_norm = 0.0;
    while (true) {
        ++step;
        this->electrochemical_evolve_one_time_step(dummy_time_step);
        solution_norm = electrochemical_solution.l2_norm();
        if (std::abs(solution_norm - old_solution_norm) < 1.0e-8) {
            break;
        } // end if
        old_solution_norm = solution_norm;
    } // end while
    if (this->verbose)
        std::cout<<step<<" iterations for finding initial electrochemical solution\n";
}

    thermal_solution = 0.0; // TODO: initialize to ambient temperature
    

    double const time_step    = input_params->get<double>("time_step"   );
    double const initial_time = input_params->get<double>("initial_time");
    double const final_time   = input_params->get<double>("final_time"  );

    std::vector<double> max_temperature;
    std::vector<double> heat_production;
    std::vector<double> voltage;
    std::vector<double> current;
    std::vector<double> time;
    std::vector<std::string> capacitor_state;
    std::vector<int> cycle;

    double current_time = initial_time;
    unsigned int step = 0;
    double data[N_DATA];

    this->thermal_setup_system(time_step);
    this->electrochemical_setup_system(time_step, GalvanostaticCharge);
    while (current_time <= final_time) {
        ++step;
        current_time += time_step;
        this->electrochemical_evolve_one_time_step(time_step);
        this->thermal_evolve_one_time_step        (time_step);
        this->process_solution(data);
        this->report_data(current_time, data);
        max_temperature.push_back(data[TEMPERATURE]);
        heat_production.push_back(data[JOULE_HEATING]);
        voltage.push_back(data[VOLTAGE]);
        current.push_back(data[CURRENT]);
        time.push_back(current_time);
        capacitor_state.push_back("charging");
        if (data[VOLTAGE] >= 2.2) {
            break;
        } //
        // TODO: ...
    } // end while

    this->electrochemical_setup_system(time_step, PotentiostaticDischarge);
    while (current_time <= final_time) {
        ++step;
        current_time += time_step;
        this->electrochemical_evolve_one_time_step(time_step);
        this->thermal_evolve_one_time_step        (time_step);
        this->process_solution(data);
        this->report_data(current_time, data);
        max_temperature.push_back(data[TEMPERATURE]);
        heat_production.push_back(data[JOULE_HEATING]);
        voltage.push_back(data[VOLTAGE]);
        current.push_back(data[CURRENT]);
        time.push_back(current_time);
        capacitor_state.push_back("discharging");
    } // end while

    output_params->put("max_temperature", to_string(max_temperature));
    output_params->put("heat_production", to_string(heat_production));
    output_params->put("voltage",         to_string(voltage));
    output_params->put("current",         to_string(current));
    output_params->put("time",            to_string(time));
    output_params->put("capacitor_state", to_string(capacitor_state));
    output_params->put("cycle",           to_string(cycle));
}                                               

template <int dim>
void 
SuperCapacitorProblem<dim>::
build_triangulation(std::shared_ptr<boost::property_tree::ptree const> database) 
{
    dealii::GridIn<dim> mesh_reader;
    std::fstream fin;
    std::string const filename = "mesh_" + std::to_string(dim) + "d.ucd";
    fin.open(filename.c_str(), std::fstream::in);
    mesh_reader.attach_triangulation(this->triangulation);
    mesh_reader.read_ucd(fin);
    fin.close();
    if (this->verbose) {
        std::cout<<"cells="<<this->triangulation.n_active_cells()<<"  "
            <<"faces="<<this->triangulation.n_active_faces()<<"  "
            <<"vertices="<<this->triangulation.n_used_vertices()<<"\n";
    }
}

template <int dim>
void SuperCapacitorProblem<dim>::set_cell_material_ids_and_face_boundary_ids(std::shared_ptr<boost::property_tree::ptree const> database) {
    double const electrode_width = 50.0e-6;
    double const separator_width = 25.0e-6;
    double const collector_width = 5.0e-6;
    double const x_min = -collector_width;
    double const x_max = 2.0 * electrode_width + separator_width + collector_width;
    double const y_min = -5.0e-6;
    double const y_max = 30.0e-6;
    double const y_top = 25.0e-6;
    double const y_bottom = 0.0e-6;
  
    dealii::types::material_id const separator_material_id         = database->get<dealii::types::material_id>("separator_material_id"        );
    dealii::types::material_id const anode_electrode_material_id   = database->get<dealii::types::material_id>("anode_electrode_material_id"  );
    dealii::types::material_id const anode_collector_material_id   = database->get<dealii::types::material_id>("anode_collector_material_id"  );
    dealii::types::material_id const cathode_electrode_material_id = database->get<dealii::types::material_id>("cathode_electrode_material_id");
    dealii::types::material_id const cathode_collector_material_id = database->get<dealii::types::material_id>("cathode_collector_material_id");

    dealii::types::boundary_id const cathode_boundary_id = database->get<dealii::types::boundary_id>("cathode_boundary_id");
    dealii::types::boundary_id const anode_boundary_id   = database->get<dealii::types::boundary_id>("anode_boundary_id"  );
    dealii::types::boundary_id const upper_boundary_id   = database->get<dealii::types::boundary_id>("upper_boundary_id"  );
    dealii::types::boundary_id const lower_boundary_id   = database->get<dealii::types::boundary_id>("lower_boundary_id"  );
    dealii::types::boundary_id const other_boundary_id   = database->get<dealii::types::boundary_id>("other_boundary_id"  );

    typename dealii::Triangulation<dim>::active_cell_iterator cell = this->triangulation.begin_active();
    typename dealii::Triangulation<dim>::active_cell_iterator end_cell = this->triangulation.end();
    for ( ; cell != end_cell; ++cell) {
        if ((cell->center()[0] > x_min)
            && (cell->center()[0] < x_min + collector_width)) {
            cell->set_material_id(anode_collector_material_id);
        } else if ((cell->center()[0] > x_min + collector_width)
            && (cell->center()[0] < x_min + collector_width + electrode_width)) {
            cell->set_material_id(anode_electrode_material_id);
        } else if ((cell->center()[0] > x_min + collector_width + electrode_width)
            && (cell->center()[0] < x_min + collector_width + electrode_width + separator_width)) {
            cell->set_material_id(separator_material_id);
        } else if ((cell->center()[0] > x_min + collector_width + electrode_width + separator_width)
            && (cell->center()[0] < x_min + collector_width + 2.0*electrode_width + separator_width)) {
            cell->set_material_id(cathode_electrode_material_id);
        } else if ((cell->center()[0] > x_min + collector_width + 2.0*electrode_width + separator_width)
            && (cell->center()[0] < x_max)) {
            cell->set_material_id(cathode_collector_material_id);
        } else {
            std::runtime_error("Error while setting material ids");
        } // end if
        for (std::size_t face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
            if (cell->face(face)->at_boundary()) {
                if ((std::abs(cell->face(face)->center()[1] - y_max) < 1.0e-10)
                    && (cell->material_id() == anode_collector_material_id)) {
                    cell->face(face)->set_boundary_indicator(anode_boundary_id);
                } else if ((std::abs(cell->face(face)->center()[1] - y_min) < 1.0e-10)
                    && (cell->material_id() == cathode_collector_material_id)) {
                    cell->face(face)->set_boundary_indicator(cathode_boundary_id);
                } else if ((std::abs(cell->face(face)->center()[1] - y_top) < 1.0e-10)
                    && ((cell->material_id() == cathode_electrode_material_id)
                        || (cell->material_id() == anode_electrode_material_id)
                        || (cell->material_id() == separator_material_id))
                    ) {
                    cell->face(face)->set_boundary_indicator(upper_boundary_id);
                } else if ((std::abs(cell->face(face)->center()[1] - y_bottom) < 1.0e-10)
                    && ((cell->material_id() == cathode_electrode_material_id)
                        || (cell->material_id() == anode_electrode_material_id)
                        || (cell->material_id() == separator_material_id))
                    ) {
                    cell->face(face)->set_boundary_indicator(lower_boundary_id);
                } else {
                    cell->face(face)->set_boundary_indicator(other_boundary_id);
                } // end if
            } // end if face at boundary
        } // end for face
    } // end for cell
}



template <int dim>
void 
SuperCapacitorProblem<dim>::
initialize_system(std::shared_ptr<boost::property_tree::ptree const> database) 
{
    // distribute degrees of freedom
    this->fe = std::make_shared<typename dealii::FESystem<dim> >
        (typename dealii::FE_Q<dim>(1), 3); // TODO: read degree p from database
    this->dof_handler.distribute_dofs(*this->fe);
    dealii::DoFRenumbering::component_wise(this->dof_handler);
    unsigned int const n_components = dealii::DoFTools::n_components(this->dof_handler);
    std::vector<dealii::types::global_dof_index> dofs_per_component(n_components);
    dealii::DoFTools::count_dofs_per_component(this->dof_handler, dofs_per_component);

    unsigned int const temperature_component      = database->get<unsigned int>("temperature_component"     );
    unsigned int const solid_potential_component  = database->get<unsigned int>("solid_potential_component" );
    unsigned int const liquid_potential_component = database->get<unsigned int>("liquid_potential_component");
    unsigned int const thermal_block              = database->get<unsigned int>("thermal_block"             );
    unsigned int const electrochemical_block      = database->get<unsigned int>("electrochemical_block"     );
    dealii::types::global_dof_index const n_thermal_dofs         = dofs_per_component[temperature_component];
    dealii::types::global_dof_index const n_electrochemical_dofs = 
        dofs_per_component[solid_potential_component] + dofs_per_component[liquid_potential_component];

    // make sparsity pattern
    unsigned int const max_couplings = this->dof_handler.max_couplings_between_dofs();
    this->sparsity_pattern.reinit(2, 2);
    this->sparsity_pattern.block(electrochemical_block, electrochemical_block).reinit(n_electrochemical_dofs, n_electrochemical_dofs, 2*max_couplings);
    this->sparsity_pattern.block(electrochemical_block, thermal_block        ).reinit(n_electrochemical_dofs, n_thermal_dofs,           max_couplings);
    this->sparsity_pattern.block(thermal_block        , electrochemical_block).reinit(n_thermal_dofs,         n_electrochemical_dofs, 2*max_couplings);
    this->sparsity_pattern.block(thermal_block        , thermal_block        ).reinit(n_thermal_dofs,         n_thermal_dofs,           max_couplings);
    this->sparsity_pattern.collect_sizes();

    dealii::DoFTools::make_sparsity_pattern(this->dof_handler, this->sparsity_pattern, this->constraint_matrix);
    this->sparsity_pattern.compress();

    // initialize matrices and vectors
    this->system_matrix.reinit(this->sparsity_pattern);
    this->system_rhs.reinit(2);
    this->system_rhs.block(electrochemical_block).reinit(n_electrochemical_dofs);
    this->system_rhs.block(thermal_block        ).reinit(n_thermal_dofs        );
    this->system_rhs.collect_sizes();

    this->solution.reinit(this->system_rhs);

    // TODO: add const ref thermal_solution and elctrochemical_solution for
    // readibility
    this->thermal_load_vector.reinit(this->system_rhs.block(thermal_block));

    if (this->verbose) {
        std::cout
            <<"total degrees of freedom : "<<this->dof_handler.n_dofs()<<"\n"
            <<"    electrochemical : "
            <<"solid_potential "<<dofs_per_component[solid_potential_component]
            <<" + liquid_potential "<<dofs_per_component[liquid_potential_component]
            <<"\n"
            <<"    thermal         : "
            <<"temperature "<<dofs_per_component[temperature_component]
            <<"\n";
    }
}

template <int dim>
void 
SuperCapacitorProblem<dim>::
reset(std::shared_ptr<boost::property_tree::ptree const> database) 
{
    std::shared_ptr<boost::property_tree::ptree> material_properties_database = 
        std::make_shared<boost::property_tree::ptree>(database->get_child("material_properties"));
    std::shared_ptr<SuperCapacitorMPValues<dim> > mp_values = std::shared_ptr<SuperCapacitorMPValues<dim> >
        (new SuperCapacitorMPValues<dim>(SuperCapacitorMPValuesParameters<dim>(material_properties_database)));

    std::shared_ptr<boost::property_tree::ptree> boundary_values_database = 
        std::make_shared<boost::property_tree::ptree>(database->get_child("boundary_values"));
    std::shared_ptr<SuperCapacitorBoundaryValues<dim> > boundary_values = std::shared_ptr<SuperCapacitorBoundaryValues<dim> >
        (new SuperCapacitorBoundaryValues<dim>(SuperCapacitorBoundaryValuesParameters<dim>(boundary_values_database)));

    // initialize operators
    this->electrochemical_operator_params = std::shared_ptr<ElectrochemicalOperatorParameters<dim> >
        (new ElectrochemicalOperatorParameters<dim>(database));
    this->electrochemical_operator_params->dof_handler       = &(this->dof_handler);
    this->electrochemical_operator_params->constraint_matrix = &(this->constraint_matrix);
    this->electrochemical_operator_params->sparsity_pattern  = &(this->sparsity_pattern.block(0, 0));
    this->electrochemical_operator_params->some_vector       = &(this->solution.block(0));
    this->electrochemical_operator_params->mp_values         = std::dynamic_pointer_cast<MPValues<dim>       const>(mp_values);
    this->electrochemical_operator_params->boundary_values   = std::dynamic_pointer_cast<BoundaryValues<dim> const>(boundary_values);
    this->electrochemical_operator = std::shared_ptr<ElectrochemicalOperator<dim> >
        (new ElectrochemicalOperator<dim>(this->electrochemical_operator_params));

    // set null space
    this->electrochemical_operator->set_null_space(database->get<unsigned int>("solid_potential_component" ), database->get<dealii::types::material_id>("material_properties.separator_material_id"        ));
    this->electrochemical_operator->set_null_space(database->get<unsigned int>("liquid_potential_component"), database->get<dealii::types::material_id>("material_properties.anode_collector_material_id"  ));
    this->electrochemical_operator->set_null_space(database->get<unsigned int>("liquid_potential_component"), database->get<dealii::types::material_id>("material_properties.cathode_collector_material_id"));
    std::vector<dealii::types::global_dof_index> const & null_space = this->electrochemical_operator->get_null_space();
    if (this->verbose) {
        std::cout<<"null space size : "<<null_space.size()<<"\n";
    }

    this->thermal_operator_params = std::shared_ptr<ThermalOperatorParameters<dim> >
        (new ThermalOperatorParameters<dim>(database));
    this->thermal_operator_params->dof_handler       = &(this->dof_handler);
    this->thermal_operator_params->constraint_matrix = &(this->constraint_matrix);
    this->thermal_operator_params->sparsity_pattern  = &(this->sparsity_pattern.block(1, 1));
    this->thermal_operator_params->some_vector       = &(this->solution.block(1));
    this->thermal_operator_params->mp_values         = std::dynamic_pointer_cast<MPValues<dim>       const>(mp_values      );
    this->thermal_operator_params->boundary_values   = std::dynamic_pointer_cast<BoundaryValues<dim> const>(boundary_values);
    this->thermal_operator = std::shared_ptr<ThermalOperator<dim> >
        (new ThermalOperator<dim>(this->thermal_operator_params));


    // initialize postprocessor
    this->post_processor_params = std::shared_ptr<SuperCapacitorPostprocessorParameters<dim> >
        (new SuperCapacitorPostprocessorParameters<dim>(database));
    this->post_processor_params->dof_handler     = &(this->dof_handler);
    this->post_processor_params->solution        = &(this->solution   );
    this->post_processor_params->mp_values       = std::dynamic_pointer_cast<MPValues<dim>       const>(mp_values      );
    this->post_processor_params->boundary_values = std::dynamic_pointer_cast<BoundaryValues<dim> const>(boundary_values);
    this->post_processor = std::shared_ptr<SuperCapacitorPostprocessor<dim> >
        (new SuperCapacitorPostprocessor<dim>(this->post_processor_params));

}

template <int dim>
void
SuperCapacitorProblem<dim>::
thermal_setup_system(double time_step)
{
    this->thermal_operator->reset(this->thermal_operator_params);
    // TODO: won't work because of dof shift
    AssertThrow((this->thermal_operator)->get_boundary_values().size() == 0,
        dealii::StandardExceptions::ExcMessage("Under construction"));

    unsigned int const thermal_block = 1; // TODO: ...
    dealii::SparseMatrix<double> & thermal_system_matrix = this->system_matrix.block(thermal_block, thermal_block);
    dealii::Vector<double> &       thermal_solution      = this->solution     .block(thermal_block);
    dealii::Vector<double> &       thermal_system_rhs    = this->system_rhs   .block(thermal_block);

    thermal_system_matrix.copy_from(this->thermal_operator->get_mass_matrix());
    thermal_system_matrix.add(time_step, this->thermal_operator->get_stiffness_matrix());

    dealii::MatrixTools::apply_boundary_values(this->thermal_operator->get_boundary_values(), thermal_system_matrix, thermal_solution, thermal_system_rhs, false);

    this->inverse_thermal_system_matrix.initialize(thermal_system_matrix);

}

template <int dim>
void SuperCapacitorProblem<dim>::
thermal_evolve_one_time_step(double const time_step) 
{   
    AssertThrow((this->thermal_operator)->get_boundary_values().size() == 0,
        dealii::StandardExceptions::ExcMessage("Under construction"));
    
    this->thermal_operator->get_mass_matrix().vmult(this->system_rhs.block(1), this->solution.block(1));
    this->thermal_load_vector = this->thermal_operator->get_load_vector();
    this->electrochemical_operator->compute_heat_source(this->solution, this->thermal_load_vector);
    this->system_rhs.block(1).add(time_step, thermal_load_vector);
    
    this->inverse_thermal_system_matrix.vmult(this->solution.block(1), this->system_rhs.block(1));
}

template <int dim>
void SuperCapacitorProblem<dim>::
electrochemical_evolve_one_time_step(double const time_step)
{ 
    std::map<dealii::types::global_dof_index, double>::const_iterator it,
        begin_it = this->electrochemical_operator->get_boundary_values().cbegin(),
        end_it = this->electrochemical_operator->get_boundary_values().cend();
    if (this->symmetric_correction) {
        for (it = begin_it ; it != end_it; ++it) {
            this->solution.block(0)[it->first] = 0.0;
        } // end for
    } // end if symmetric correction
        (this->electrochemical_operator)->get_mass_matrix().vmult(this->system_rhs.block(0), this->solution.block(0));
    for (it = begin_it ; it != end_it; ++it) {
        this->system_rhs.block(0)[it->first] = 0.0;
        this->solution.block(0)[it->first] = it->second;
    } // end for
    this->system_rhs.block(0).add(time_step, this->electrochemical_operator->get_load_vector());

    begin_it = this->rhs_set.cbegin();
    end_it = this->rhs_set.cend();
    for (it = begin_it ; it != end_it; ++it) {
        this->system_rhs.block(0)[it->first] = it->second;
    } // end for
    if (this->symmetric_correction) {
        begin_it = this->rhs_add.cbegin();
        end_it = this->rhs_add.cend();
        for (it = begin_it ; it != end_it; ++it) {
            this->system_rhs.block(0)[it->first] += it->second;
        } // end for
    } // end if symmetric correction

    this->inverse_electrochemical_system_matrix.vmult(this->solution.block(0), this->system_rhs.block(0));
}

template <int dim>
void
SuperCapacitorProblem<dim>::
electrochemical_setup_system(double const time_step, CapacitorState const capacitor_state) 
{   
    this->electrochemical_operator_params->capacitor_state = capacitor_state;
    this->electrochemical_operator->reset(this->electrochemical_operator_params);
        
    this->system_matrix.block(0, 0).copy_from(this->electrochemical_operator->get_mass_matrix());
    this->system_matrix.block(0, 0).add(time_step, (this->electrochemical_operator)->get_stiffness_matrix());
    
    std::vector<dealii::types::global_dof_index> const & null_space = this->electrochemical_operator->get_null_space();
    this->system_matrix.block(0, 0).set(null_space, dealii::FullMatrix<double>(dealii::IdentityMatrix(null_space.size())));
    
    this->system_rhs.block(0) = 0.0;
    dealii::MatrixTools::apply_boundary_values(this->electrochemical_operator->get_boundary_values(), this->system_matrix.block(0, 0), this->solution.block(0), this->system_rhs.block(0), this->symmetric_correction);

    this->inverse_electrochemical_system_matrix.initialize(this->system_matrix.block(0, 0));
        
    this->rhs_set.clear();
    this->rhs_add.clear();
    std::map<dealii::types::global_dof_index, double>::const_iterator it = this->electrochemical_operator->get_boundary_values().cbegin();
    std::map<dealii::types::global_dof_index, double>::const_iterator end_it = this->electrochemical_operator->get_boundary_values().cend();
    for ( ; it != end_it; ++it) {
        this->rhs_set[it->first] = this->system_rhs.block(0)[it->first];
        this->system_rhs.block(0)[it->first] = 0.0;
    } // end for
    AssertThrow(rhs_set.size() == this->electrochemical_operator->get_boundary_values().size(),
        dealii::StandardExceptions::ExcDimensionMismatch(rhs_set.size(), this->electrochemical_operator->get_boundary_values().size()));
    if ((this->symmetric_correction)
        && (this->system_rhs.block(0).l2_norm() != 0.0)) {
        dealii::types::global_dof_index n = this->system_rhs.block(0).size();
        for (dealii::types::global_dof_index i = 0; i < n; ++i) {
            if (this->system_rhs.block(0)[i] != 0.0) {
                this->rhs_add[i] = this->system_rhs.block(0)[i];
            } // end if entry is non zero
        } // end for all entries
    } // end if symetric correction
//    if (this->verbose) {
//        std::cout<<"rhs set size is "<<rhs_set.size()<<"\n"
//            "rhs add size is "<<rhs_add.size()<<"\n";
//    }
}

template <int dim>
void
SuperCapacitorProblem<dim>::
process_solution(double * data) 
{   
    this->post_processor->reset(this->post_processor_params);
    this->post_processor->get("voltage",         data[VOLTAGE]      );
    this->post_processor->get("current",         data[CURRENT]      );
    this->post_processor->get("joule_heating",   data[JOULE_HEATING]);
    this->post_processor->get("max_temperature", data[TEMPERATURE]  );
    this->post_processor->get("volume",          data[VOLUME]       );
    this->post_processor->get("mass",            data[MASS]         );
    this->post_processor->get("surface_area",    data[SURFACE_AREA] );
}          

} // end namespace cap
