#ifndef CACHE_SUPER_CAPACITOR_H
#define CACHE_SUPER_CAPACITOR_H

#include <cache/utils.h>
#include <cache/mp_values.h>
#include <cache/boundary_values.h>
#include <cache/thermal_operator.h>
#include <cache/thermal_operator.h>
#include <cache/electrochemical_operator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace cache {

class SuperCapacitorProblemParameters {
public:
    SuperCapacitorProblemParameters(std::shared_ptr<boost::property_tree::ptree const> d)
        : database(d)
        { }
    std::shared_ptr<boost::property_tree::ptree const> database;
};


template <int dim>
class SuperCapacitorProblem {
public:
    SuperCapacitorProblem(SuperCapacitorProblemParameters problem_params);
    void run(std::shared_ptr<boost::property_tree::ptree const> input_params,
             std::shared_ptr<boost::property_tree::ptree>       output_params);
                                        
private:
    void build_triangulation();
    void set_cell_material_ids_and_face_boundary_ids();
    void initialize_system();
    void reset(std::shared_ptr<boost::property_tree::ptree const> params);

    void electrochemical_setup_system(double const time_step);
    void electrochemical_evolve_one_time_step(double const time_step);
    void thermal_setup_system(double const time_step);
    void thermal_evolve_one_time_step(double const time_step);

    void process_solution();

    typename dealii::Triangulation<dim> triangulation;

    unsigned int const n_components;
    unsigned int const solid_potential_component;
    unsigned int const liquid_potential_component;
    unsigned int const temperature_component;
    unsigned int const degree_p;
    typename dealii::FESystem<dim> fe;
    typename dealii::DoFHandler<dim> dof_handler;
    dealii::ConstraintMatrix constraint_matrix;
    dealii::BlockSparsityPattern sparsity_pattern;
    dealii::BlockSparseMatrix<double> system_matrix;
    dealii::BlockVector<double> system_rhs;
    dealii::BlockVector<double> solution;

    dealii::SparseDirectUMFPACK inverse_electrochemical_system_matrix;
    dealii::SparseDirectUMFPACK inverse_thermal_system_matrix;
    dealii::Vector<double> thermal_load_vector;
    
    std::shared_ptr<cache::ElectrochemicalOperatorParameters<dim> > electrochemical_operator_params;
    std::shared_ptr<cache::ElectrochemicalOperator<dim> >           electrochemical_operator;
    std::shared_ptr<cache::ThermalOperatorParameters<dim> >         thermal_operator_params;
    std::shared_ptr<cache::ThermalOperator<dim> >                   thermal_operator;

    std::shared_ptr<SuperCapacitorMPValues<dim> >       mp_values;
    std::shared_ptr<SuperCapacitorBoundaryValues<dim> > bp_values;

    bool const symmetric_correction;
    std::map<dealii::types::global_dof_index, double> rhs_set;
    std::map<dealii::types::global_dof_index, double> rhs_add;
// TODO: ...
    std::shared_ptr<boost::property_tree::ptree const> database;
};

} // end namespace

#endif // CACHE_SUPER_CAPACITOR_H
