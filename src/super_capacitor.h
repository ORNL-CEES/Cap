#ifndef CAP_SUPER_CAPACITOR_H
#define CAP_SUPER_CAPACITOR_H

#include <cap/utils.h>
#include <cap/mp_values.h>
#include <cap/boundary_values.h>
#include <cap/thermal_operator.h>
#include <cap/electrochemical_operator.h>
#include <cap/post_processor.h>
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

namespace cap {

template <int dim>
class SuperCapacitorProblem {
public:
    SuperCapacitorProblem(std::shared_ptr<boost::property_tree::ptree const> database);
    void run(std::shared_ptr<boost::property_tree::ptree const> input_params,
             std::shared_ptr<boost::property_tree::ptree>       output_params);
                                        
private:
    void run_constant_current_charge_constant_voltage_discharge(
        std::shared_ptr<boost::property_tree::ptree const> input_params,
        std::shared_ptr<boost::property_tree::ptree>       output_params);
    void run_constant_current_cycling(
        std::shared_ptr<boost::property_tree::ptree const> input_params,
        std::shared_ptr<boost::property_tree::ptree>       output_params);
    void run_constant_power_cycling(
        std::shared_ptr<boost::property_tree::ptree const> input_params,
        std::shared_ptr<boost::property_tree::ptree>       output_params);
    void build_triangulation(std::shared_ptr<boost::property_tree::ptree const> database);
    void set_cell_material_ids_and_face_boundary_ids(std::shared_ptr<boost::property_tree::ptree const> database);
    void initialize_system(std::shared_ptr<boost::property_tree::ptree const> database);

    void reset(std::shared_ptr<boost::property_tree::ptree const> params);

    void electrochemical_setup_system(double const time_step, CapacitorState const capacitor_state);
    void electrochemical_evolve_one_time_step(double const time_step);
    void thermal_setup_system(double const time_step);
    void thermal_evolve_one_time_step(double const time_step);

    enum OutputData { TEMPERATURE, VOLTAGE, CURRENT, JOULE_HEATING, SURFACE_AREA, VOLUME, MASS, N_DATA};
    void process_solution(double * data);
    void report_data (double const time, double const * data) const;

    typename dealii::Triangulation<dim> triangulation;

    std::shared_ptr<typename dealii::FESystem<dim> > fe; // TODO: would be nice to get rid of this guy
    typename dealii::DoFHandler<dim>    dof_handler;
    dealii::ConstraintMatrix            constraint_matrix;
    dealii::BlockSparsityPattern        sparsity_pattern;
    dealii::BlockSparseMatrix<double>   system_matrix;
    dealii::BlockVector<double>         system_rhs;
    dealii::BlockVector<double>         solution;

    dealii::SparseDirectUMFPACK inverse_electrochemical_system_matrix;
    dealii::SparseDirectUMFPACK inverse_thermal_system_matrix;
    dealii::Vector<double> thermal_load_vector;
    
    std::shared_ptr<ElectrochemicalOperatorParameters<dim> > electrochemical_operator_params;
    std::shared_ptr<ElectrochemicalOperator<dim> >           electrochemical_operator;
    std::shared_ptr<ThermalOperatorParameters<dim> >         thermal_operator_params;
    std::shared_ptr<ThermalOperator<dim> >                   thermal_operator;

    std::shared_ptr<SuperCapacitorPostprocessorParameters<dim> > post_processor_params;
    std::shared_ptr<SuperCapacitorPostprocessor<dim> >           post_processor;

    bool const verbose;

    bool const symmetric_correction;
    std::map<dealii::types::global_dof_index, double> rhs_set;
    std::map<dealii::types::global_dof_index, double> rhs_add;
};

} // end namespace cap

#endif // CAP_SUPER_CAPACITOR_H
