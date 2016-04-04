/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_DEAL_II_SUPERCAPACITOR_H
#define CAP_DEAL_II_SUPERCAPACITOR_H

#include <cap/energy_storage_device.h>
#include <cap/geometry.h>
#include <cap/electrochemical_physics.h>
#include <cap/post_processor.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/block_vector.h>
#include <memory>
#include <iostream>

namespace cap
{

template <int dim>
class SuperCapacitorInspector : public EnergyStorageDeviceInspector
{
public:
  SuperCapacitorInspector() = default;
  void inspect(EnergyStorageDevice *device);
};

template <int dim>
class SuperCapacitor : public EnergyStorageDevice
{
public:
  SuperCapacitor(boost::mpi::communicator const &comm,
                 boost::property_tree::ptree const &ptree);

  void inspect(EnergyStorageDeviceInspector *inspector) override;

  void get_voltage(double &voltage) const override;

  void get_current(double &current) const override;

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
   * Granting access to the post-processor for the inspector.
   */
  std::shared_ptr<Postprocessor<dim>> get_post_processor() const;

  /**
   * Provides a copy of the property tree used to build the supercapacitor to
   * the inspector.
   */
  boost::property_tree::ptree const *get_property_tree() const;

private:
  /**
   * Helper function to advance time by @p time_step second.
   */
  void evolve_one_time_step(double const time_step,
                            SuperCapacitorState supercapacitor_state,
                            bool rebuild);

  /**
   * Maximum number of iterations of the Krylov solver in
   * evolve_one_time_step().
   */
  unsigned int max_iter;
  /**
   * Verbosity level of the Krylov solver in evolve_one_time_step().
   */
  unsigned int verbose_lvl;
  /**
   * Absolute tolerance of the Krylov solver in evolve_one_time_step(). The
   * tolerance used by the Krylov solver is the maximum of the relative and the
   * absolute tolerance.
   */
  double abs_tolerance;
  /**
   * Relative tolerance of the Krylov solver in evolve_one_time_step(), i.e. the
   * tolerance is @p rel_toleracne \f$ \times ||b||_{2}\f$. The tolerance used
   * by the Krylov solver is the maximum of the relative and the absolute
   * tolerance.
   */
  double rel_tolerance;
  /**
   * Area of the cathode.
   */
  double surface_area;

  std::shared_ptr<SuperCapacitorGeometry<dim>> geometry;
  std::shared_ptr<dealii::FESystem<dim>> fe;
  std::shared_ptr<dealii::DoFHandler<dim>> dof_handler;
  std::shared_ptr<dealii::BlockVector<double>> solution;

  std::shared_ptr<ElectrochemicalPhysicsParameters<dim>>
      electrochemical_physics_params;
  std::shared_ptr<ElectrochemicalPhysics<dim>> electrochemical_physics;
  std::shared_ptr<SuperCapacitorPostprocessorParameters<dim>>
      post_processor_params;
  std::shared_ptr<SuperCapacitorPostprocessor<dim>> post_processor;
  boost::property_tree::ptree const _ptree;

  template <int dimension>
  friend class SuperCapacitorInspector;
};
}

#endif
