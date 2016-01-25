#ifndef CAP_SUPER_CAPACITOR_H
#define CAP_SUPER_CAPACITOR_H

#include <cap/energy_storage_device.h>
#include <cap/geometry.h>
#include <cap/thermal_operator.h>
#include <cap/electrochemical_operator.h>
#include <cap/post_processor.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <memory>
#include <iostream>

namespace cap
{

template <int dim>
class SuperCapacitor : public EnergyStorageDevice
{
public:
  SuperCapacitor(boost::mpi::communicator const &comm,
                 boost::property_tree::ptree const &ptree);
  void inspect(EnergyStorageDeviceInspector *inspector) override;
  void print_data(std::ostream &os) const override;
  void get_voltage(double &voltage) const override;
  void get_current(double &current) const override;
  void reset_voltage(double const voltage) override;
  void reset_current(double const current) override;
  void
  evolve_one_time_step_constant_current(double const time_step,
                                        double const constant_current) override;
  void
  evolve_one_time_step_constant_voltage(double const time_step,
                                        double const constant_voltage) override;
  void
  evolve_one_time_step_constant_power(double const time_step,
                                      double const constant_power) override;
  void evolve_one_time_step_constant_load(double const time_step,
                                          double const constant_load) override;

  void evolve_one_time_step(double const time_step);

private:
  std::shared_ptr<dealii::FESystem<dim>>
      fe; // TODO: would be nice to get rid of this guy
  std::shared_ptr<dealii::DoFHandler<dim>> dof_handler;
  std::shared_ptr<dealii::ConstraintMatrix> constraint_matrix;
  std::shared_ptr<dealii::BlockSparsityPattern> sparsity_pattern;
  std::shared_ptr<dealii::BlockSparseMatrix<double>> system_matrix;
  std::shared_ptr<dealii::BlockVector<double>> system_rhs;
  std::shared_ptr<dealii::BlockVector<double>> solution;

  //    dealii::SparseDirectUMFPACK inverse_electrochemical_system_matrix;
  //    dealii::SparseDirectUMFPACK inverse_thermal_system_matrix;
  //    dealii::Vector<double> thermal_load_vector;
  //
  std::shared_ptr<SuperCapacitorGeometry<dim>> geometry;
  std::shared_ptr<ElectrochemicalOperatorParameters<dim>>
      electrochemical_operator_params;
  std::shared_ptr<ElectrochemicalOperator<dim>> electrochemical_operator;
  std::shared_ptr<ThermalOperatorParameters<dim>> thermal_operator_params;
  std::shared_ptr<ThermalOperator<dim>> thermal_operator;
  std::shared_ptr<SuperCapacitorPostprocessorParameters<dim>>
      post_processor_params;
  std::shared_ptr<SuperCapacitorPostprocessor<dim>> post_processor;
};

} // end namespace cap

#endif // CAP_SUPER_CAPACITOR_H
