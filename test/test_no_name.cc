#define BOOST_TEST_MODULE NoName
#define BOOST_TEST_MAIN
#include <cap/geometry.h>
#include <cap/energy_storage_device.h>
#include <cap/thermal_operator.h>
#include <cap/electrochemical_operator.h>
#include <cap/post_processor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/format.hpp>
#include <memory>
#include <tuple>
#include <fstream>



namespace cap {

template <int dim>
class NoName : public EnergyStorageDevice
{
public:
    NoName(std::shared_ptr<Parameters const> parameters);
    void print_data(std::ostream & os) const;
    void reset_voltage(double const voltage);
    void reset_current(double const current);
    void evolve_one_time_step_constant_current(double const time_step, double const constant_current);
    void evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage);
    void evolve_one_time_step_constant_power  (double const time_step, double const constant_power  );

    void evolve_one_time_step(double const time_step);
private:
    std::shared_ptr<typename dealii::FESystem<dim>    > fe; // TODO: would be nice to get rid of this guy
    std::shared_ptr<typename dealii::DoFHandler<dim>  > dof_handler;
    std::shared_ptr<dealii::ConstraintMatrix          > constraint_matrix;
    std::shared_ptr<dealii::BlockSparsityPattern      > sparsity_pattern;
    std::shared_ptr<dealii::BlockSparseMatrix<double> > system_matrix;
    std::shared_ptr<dealii::BlockVector<double>       > system_rhs;
    std::shared_ptr<dealii::BlockVector<double>       > solution;

//    dealii::SparseDirectUMFPACK inverse_electrochemical_system_matrix;
//    dealii::SparseDirectUMFPACK inverse_thermal_system_matrix;
//    dealii::Vector<double> thermal_load_vector;
//                                                               
    std::shared_ptr<cap::SuperCapacitorGeometry<dim>           > geometry;

    std::shared_ptr<ElectrochemicalOperatorParameters<dim>     > electrochemical_operator_params;
    std::shared_ptr<ElectrochemicalOperator<dim>               > electrochemical_operator;
    std::shared_ptr<ThermalOperatorParameters<dim>             > thermal_operator_params;
    std::shared_ptr<ThermalOperator<dim>                       > thermal_operator;

    std::shared_ptr<SuperCapacitorPostprocessorParameters<dim> > post_processor_params;
    std::shared_ptr<SuperCapacitorPostprocessor<dim>           > post_processor;
};



template <int dim>
NoName<dim>::
NoName(std::shared_ptr<Parameters const> parameters)
{
    // get database
    std::shared_ptr<boost::property_tree::ptree const> database = parameters->database;

    // build triangulation
    this->geometry = std::make_shared<cap::SuperCapacitorGeometry<dim> >(database);
    std::shared_ptr<dealii::Triangulation<dim> const> triangulation = (*this->geometry).get_triangulation();
    
    // distribute degrees of freedom
    this->fe          = std::make_shared<dealii::FESystem<dim>   >(dealii::FE_Q<dim>(1), 3);
    this->dof_handler = std::make_shared<dealii::DoFHandler<dim> >(*triangulation);
    (*this->dof_handler).distribute_dofs(*this->fe);
    dealii::DoFRenumbering::component_wise(*this->dof_handler);
    unsigned int const n_components = dealii::DoFTools::n_components(*this->dof_handler);
    std::vector<dealii::types::global_dof_index> dofs_per_component(n_components);
    dealii::DoFTools::count_dofs_per_component(*this->dof_handler, dofs_per_component);

    unsigned int const temperature_component      = database->get<unsigned int>("temperature_component"     );
    unsigned int const solid_potential_component  = database->get<unsigned int>("solid_potential_component" );
    unsigned int const liquid_potential_component = database->get<unsigned int>("liquid_potential_component");
//    unsigned int const n_components               = database->get<unsigned int>("n_components"              );
    unsigned int const thermal_block              = database->get<unsigned int>("thermal_block"             );
    unsigned int const electrochemical_block      = database->get<unsigned int>("electrochemical_block"     );
    unsigned int const n_blocks                   = database->get<unsigned int>("n_blocks"                  );
    dealii::types::global_dof_index const n_thermal_dofs         =
        dofs_per_component[temperature_component];
    dealii::types::global_dof_index const n_electrochemical_dofs =
        dofs_per_component[solid_potential_component] + dofs_per_component[liquid_potential_component];

    // make sparsity pattern
    this->sparsity_pattern  = std::make_shared<dealii::BlockSparsityPattern>();
    this->constraint_matrix = std::make_shared<dealii::ConstraintMatrix    >();
    unsigned int const max_couplings = (*this->dof_handler).max_couplings_between_dofs();
    (*this->sparsity_pattern).reinit(n_blocks, n_blocks);
    (*this->sparsity_pattern).block(electrochemical_block, electrochemical_block).reinit(n_electrochemical_dofs, n_electrochemical_dofs, 2*max_couplings);
    (*this->sparsity_pattern).block(electrochemical_block, thermal_block        ).reinit(n_electrochemical_dofs, n_thermal_dofs,           max_couplings);
    (*this->sparsity_pattern).block(thermal_block        , electrochemical_block).reinit(n_thermal_dofs,         n_electrochemical_dofs, 2*max_couplings);
    (*this->sparsity_pattern).block(thermal_block        , thermal_block        ).reinit(n_thermal_dofs,         n_thermal_dofs,           max_couplings);
    (*this->sparsity_pattern).collect_sizes();
    dealii::DoFTools::make_sparsity_pattern(*this->dof_handler, *this->sparsity_pattern, *this->constraint_matrix);
    (*this->sparsity_pattern).compress();

    // initialize matrices and vectors
    this->system_matrix = std::make_shared<dealii::BlockSparseMatrix<double> >();
    this->system_rhs    = std::make_shared<dealii::BlockVector<double>       >();
    this->solution      = std::make_shared<dealii::BlockVector<double>       >();
    (*this->system_matrix).reinit(*this->sparsity_pattern);
    (*this->system_rhs).reinit(n_blocks);
    (*this->system_rhs).block(electrochemical_block).reinit(n_electrochemical_dofs);
    (*this->system_rhs).block(thermal_block        ).reinit(n_thermal_dofs        );
    (*this->system_rhs).collect_sizes();
    (*this->solution).reinit(*this->system_rhs);



    std::shared_ptr<boost::property_tree::ptree> material_properties_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child("material_properties"));
    std::shared_ptr<SuperCapacitorMPValues<dim> > mp_values =
        std::make_shared<SuperCapacitorMPValues<dim> >(SuperCapacitorMPValuesParameters<dim>(material_properties_database));

    std::shared_ptr<boost::property_tree::ptree> boundary_values_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child("boundary_values")    );
    std::shared_ptr<SuperCapacitorBoundaryValues<dim> > boundary_values =
        std::make_shared<SuperCapacitorBoundaryValues<dim> >(SuperCapacitorBoundaryValuesParameters<dim>(boundary_values_database));

    // initialize electrochemical operator
    this->electrochemical_operator_params =
        std::make_shared<ElectrochemicalOperatorParameters<dim> >(database);
    this->electrochemical_operator_params->dof_handler       = this->dof_handler;
    this->electrochemical_operator_params->constraint_matrix = this->constraint_matrix;
    this->electrochemical_operator_params->sparsity_pattern  = &((*this->sparsity_pattern).block(electrochemical_block, electrochemical_block));
    this->electrochemical_operator_params->some_vector       = &((*this->solution).block(electrochemical_block));
    this->electrochemical_operator_params->mp_values         = std::dynamic_pointer_cast<MPValues<dim>       const>(mp_values);
    this->electrochemical_operator_params->boundary_values   = std::dynamic_pointer_cast<BoundaryValues<dim> const>(boundary_values);
    this->electrochemical_operator =
        std::make_shared<ElectrochemicalOperator<dim> >(this->electrochemical_operator_params);

    // set null space
    dealii::types::material_id separator_material_id         = database->get<dealii::types::material_id>("material_properties.separator_material_id"        );
    dealii::types::material_id anode_collector_material_id   = database->get<dealii::types::material_id>("material_properties.anode_collector_material_id"  );
    dealii::types::material_id cathode_collector_material_id = database->get<dealii::types::material_id>("material_properties.cathode_collector_material_id");
    this->electrochemical_operator->set_null_space(solid_potential_component , separator_material_id        );
    this->electrochemical_operator->set_null_space(liquid_potential_component, anode_collector_material_id  );
    this->electrochemical_operator->set_null_space(liquid_potential_component, cathode_collector_material_id);

    // initialize thermal operator
    this->thermal_operator_params =
        std::make_shared<ThermalOperatorParameters<dim> >(database);
    this->thermal_operator_params->dof_handler       = this->dof_handler;
    this->thermal_operator_params->constraint_matrix = this->constraint_matrix;
    this->thermal_operator_params->sparsity_pattern  = &((*this->sparsity_pattern).block(thermal_block, thermal_block));
    this->thermal_operator_params->some_vector       = &((*this->solution).block(thermal_block));
    this->thermal_operator_params->mp_values         = std::dynamic_pointer_cast<MPValues<dim>       const>(mp_values      );
    this->thermal_operator_params->boundary_values   = std::dynamic_pointer_cast<BoundaryValues<dim> const>(boundary_values);
    this->thermal_operator =
        std::make_shared<ThermalOperator<dim> >(this->thermal_operator_params);

    // initialize postprocessor
    this->post_processor_params =
        std::make_shared<SuperCapacitorPostprocessorParameters<dim> >(database);
    this->post_processor_params->dof_handler     = this->dof_handler;
    this->post_processor_params->solution        = this->solution.get(); // TODO:
    this->post_processor_params->mp_values       = std::dynamic_pointer_cast<MPValues<dim>       const>(mp_values      );
    this->post_processor_params->boundary_values = std::dynamic_pointer_cast<BoundaryValues<dim> const>(boundary_values);
    this->post_processor =
        std::make_shared<SuperCapacitorPostprocessor<dim> >(this->post_processor_params);
}



template <int dim>
void
NoName<dim>::
evolve_one_time_step_constant_voltage(double const time_step, double const constant_voltage)
{
    this->electrochemical_operator_params->capacitor_state = CustomConstantVoltage;
    this->electrochemical_operator_params->custom_constant_voltage = constant_voltage;
    this->evolve_one_time_step(time_step);
}



template <int dim>
void
NoName<dim>::
evolve_one_time_step_constant_current(double const time_step, double const constant_current)
{
    this->electrochemical_operator_params->capacitor_state = CustomConstantCurrent;
    this->electrochemical_operator_params->custom_constant_current = constant_current;
    this->evolve_one_time_step(time_step);
}

template <int dim>
void
NoName<dim>::
evolve_one_time_step(double const time_step)
{
    bool const symmetric_correction = true;
    this->electrochemical_operator->reset(this->electrochemical_operator_params);
        
    unsigned int const electrochemical_block = 1; // TODO: !!!!!!!!
    dealii::SparseMatrix<double> & system_matrix_electrochemical_block = (this->system_matrix)->block(electrochemical_block, electrochemical_block);
    dealii::Vector      <double> & system_rhs_electrochemical_block    = (this->system_rhs   )->block(electrochemical_block);  
    dealii::Vector      <double> & solution_electrochemical_block      = (this->solution     )->block(electrochemical_block);  
    system_matrix_electrochemical_block.copy_from((this->electrochemical_operator)->get_mass_matrix());
    system_matrix_electrochemical_block.add(time_step, (this->electrochemical_operator)->get_stiffness_matrix());
    
    std::vector<dealii::types::global_dof_index> const & null_space = (this->electrochemical_operator)->get_null_space();
    system_matrix_electrochemical_block.set(null_space, dealii::FullMatrix<double>(dealii::IdentityMatrix(null_space.size())));
    
    system_rhs_electrochemical_block = 0.0;
    dealii::MatrixTools::apply_boundary_values(this->electrochemical_operator->get_boundary_values(), system_matrix_electrochemical_block, solution_electrochemical_block, system_rhs_electrochemical_block, symmetric_correction);

    dealii::SparseDirectUMFPACK inverse_electrochemical_system_matrix;
    inverse_electrochemical_system_matrix.initialize(system_matrix_electrochemical_block);
      
    std::map<dealii::types::global_dof_index, double> rhs_set;
    std::map<dealii::types::global_dof_index, double> rhs_add;
    rhs_set.clear();
    rhs_add.clear();
    std::map<dealii::types::global_dof_index, double>::const_iterator it;
    std::map<dealii::types::global_dof_index, double>::const_iterator begin_it = (this->electrochemical_operator)->get_boundary_values().cbegin();
    std::map<dealii::types::global_dof_index, double>::const_iterator end_it   = (this->electrochemical_operator)->get_boundary_values().cend();
    for (it = begin_it ; it != end_it; ++it) {
        rhs_set[it->first] = system_rhs_electrochemical_block[it->first];
        system_rhs_electrochemical_block[it->first] = 0.0;
    } // end for
    AssertThrow(rhs_set.size() == this->electrochemical_operator->get_boundary_values().size(),
        dealii::StandardExceptions::ExcDimensionMismatch(rhs_set.size(), this->electrochemical_operator->get_boundary_values().size()));
    if ((symmetric_correction)
        && (system_rhs_electrochemical_block.l2_norm() != 0.0)) {
        dealii::types::global_dof_index n = system_rhs_electrochemical_block.size();
        for (dealii::types::global_dof_index i = 0; i < n; ++i) {
            if (system_rhs_electrochemical_block[i] != 0.0) {
                rhs_add[i] = system_rhs_electrochemical_block[i];
            } // end if entry is non zero
        } // end for all entries
    } // end if symetric correction



    begin_it = (this->electrochemical_operator)->get_boundary_values().cbegin(),
    end_it   = (this->electrochemical_operator)->get_boundary_values().cend();
    if (symmetric_correction) {
        for (it = begin_it ; it != end_it; ++it) {
            solution_electrochemical_block[it->first] = 0.0;
        } // end for
    } // end if symmetric correction
        (this->electrochemical_operator)->get_mass_matrix().vmult(system_rhs_electrochemical_block, solution_electrochemical_block);
    for (it = begin_it ; it != end_it; ++it) {
        system_rhs_electrochemical_block[it->first] = 0.0;
        solution_electrochemical_block[it->first] = it->second;
    } // end for
    system_rhs_electrochemical_block.add(time_step, this->electrochemical_operator->get_load_vector());

    begin_it = rhs_set.cbegin();
    end_it   = rhs_set.cend();
    for (it = begin_it ; it != end_it; ++it) {
        system_rhs_electrochemical_block[it->first] = it->second;
    } // end for
    if (symmetric_correction) {
        begin_it = rhs_add.cbegin();
        end_it   = rhs_add.cend();
        for (it = begin_it ; it != end_it; ++it) {
            system_rhs_electrochemical_block[it->first] += it->second;
        } // end for
    } // end if symmetric correction

    inverse_electrochemical_system_matrix.vmult(solution_electrochemical_block, system_rhs_electrochemical_block);
}



template <int dim>
void
NoName<dim>::
evolve_one_time_step_constant_power(double const time_step, double const constant_power)
{
    std::ignore = time_step;
    std::ignore = constant_power;
    throw std::runtime_error("not implemented");
}



template <int dim>
void
NoName<dim>::
print_data(std::ostream & os) const
{
    double current;
    double voltage;
    this->post_processor->reset(this->post_processor_params);
    this->post_processor->get("voltage", voltage);
    this->post_processor->get("current", current);
    os<<boost::format("  %10.5f  %10.7f  \n")
        % current
        % voltage
        ;
}



template <int dim>
void
NoName<dim>::
reset_voltage(double const voltage)
{
    std::ignore = voltage;
    throw std::runtime_error("not implemented");
}



template <int dim>
void
NoName<dim>::
reset_current(double const current)
{
    std::ignore = current;
    throw std::runtime_error("not implemented");
}




} // end namespace cap

BOOST_AUTO_TEST_CASE( test_no_name )
{
    std::shared_ptr<boost::property_tree::ptree> database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_no_name", *database);

    cap::NoName<2> no_name(std::make_shared<cap::Parameters const>(database));

    double const charge_time     = database->get<double>("charge_time"    );
    double const relaxation_time = database->get<double>("relaxation_time");
    double const time_step       = database->get<double>("time_step"      );
    double const charge_current  = database->get<double>("charge_current" );
    std::fstream fout;
    fout.open("no_name_data", std::fstream::out);
    std::ostream & os = fout;
    for (double time = 0.0 ; time < charge_time; time += time_step)
    {
        no_name.evolve_one_time_step_constant_current(time_step, charge_current);
        os<<boost::format("%5.1f  ") % time;
        no_name.print_data(os);
    }
    for (double time = charge_time ; time < charge_time+relaxation_time; time += time_step)
    {
        no_name.evolve_one_time_step_constant_current(time_step, 0.0);
        os<<boost::format("%5.1f  ") % time;
        no_name.print_data(os);
    }
}
