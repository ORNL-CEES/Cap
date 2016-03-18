/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

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
  [[deprecated]] SuperCapacitor(boost::mpi::communicator const &comm,
                                boost::property_tree::ptree const &ptree);

  void inspect(EnergyStorageDeviceInspector *inspector) override;

  void print_data(std::ostream &os) const override;

  void get_voltage(double &voltage) const override;

  void get_current(double &current) const override;

  void reset_voltage(double const voltage) override;

  void reset_current(double const current) override;

  void evolve_one_time_step_constant_current(double const time_step,
                                             double const current) override;

  void evolve_one_time_step_constant_voltage(double const time_step,
                                             double const voltage) override;

  void evolve_one_time_step_constant_power(double const time_step,
                                           double const power) override;

  void evolve_one_time_step_constant_load(double const time_step,
                                          double const load) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_current(double const time_step,
                                           double const current) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_voltage(double const time_step,
                                           double const voltage) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_power(double const time_step,
                                         double const power) override;

  /**
   * This function is not implemented and throws an exception.
   */
  void evolve_one_time_step_linear_load(double const time_step,
                                        double const load) override;

  /**
   * Helper function to advance time by @p time_step second.
   */
  // TODO: move to private.
  void evolve_one_time_step(double const time_step);

private:
  /**
   * Maximum number of iterations of the Krylov solver in
   * evolve_one_time_step().
   */
  unsigned int max_iter;
  /**
   * Relative tolerance of the Krylov solver in evolve_one_time_step(), i.e. the
   * tolerance is @p rel_toleracne \f$ \times ||b||_{2}\f$. The tolerance used
   * by the Krylov solver is the maximum of the relative and the absolute
   * tolerance.
   */
  double rel_tolerance;
  /**
   * Absolute tolerance of the Krylov solver in evolve_one_time_step(). The
   * tolerance used by the Krylov solver is the maximum of the relative and the
   * absolute tolerance.
   */
  double abs_tolerance;
  /**
   * Verbosity level of the Krylov solver in evolve_one_time_step().
   */
  unsigned int verbose_lvl;

  std::shared_ptr<SuperCapacitorGeometry<dim>> geometry;
  // TODO: would be nice to get rid of this guy
  std::shared_ptr<dealii::FESystem<dim>> fe;
  std::shared_ptr<dealii::DoFHandler<dim>> dof_handler;
  std::shared_ptr<dealii::ConstraintMatrix> constraint_matrix;
  std::shared_ptr<dealii::BlockSparsityPattern> sparsity_pattern;
  std::shared_ptr<dealii::BlockSparseMatrix<double>> system_matrix;
  std::shared_ptr<dealii::BlockVector<double>> system_rhs;
  std::shared_ptr<dealii::BlockVector<double>> solution;

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
