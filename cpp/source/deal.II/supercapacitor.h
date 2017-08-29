/* Copyright (c) 2016 - 2017, the Cap authors.
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
#include <cap/timer.h>
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

  /**
   * Output the mesh with the subdomain IDs associated to each cells. It can
   * also be used to output some quantities of interest.
   */
  void inspect(EnergyStorageDevice *device) override;
};

template <int dim>
class SuperCapacitor : public EnergyStorageDevice
{
public:
  SuperCapacitor(boost::property_tree::ptree const &ptree,
                 boost::mpi::communicator const &comm);

  ~SuperCapacitor();

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
   * Return the underlying Geometry.
   */
  std::shared_ptr<Geometry<dim>> get_geometry() const;

  /**
   * Granting access to the post-processor parameters for the inspector.
   */
  std::shared_ptr<PostprocessorParameters<dim>>
  get_post_processor_parameters() const;

  /**
   * Granting access to the post-processor for the inspector.
   */
  std::shared_ptr<Postprocessor<dim>> get_post_processor() const;

  /**
   * Provides a copy of the property tree used to build the supercapacitor to
   * the inspector.
   */
  boost::property_tree::ptree const *get_property_tree() const;

  /**
   * Save the current state of energy device in a compressed file.
   */
  void save(const std::string &filename) const override;

  /**
   * Load an energy device from a state saved in a compressed file.
   */
  void load(const std::string &filename) override;

private:
  /**
   * Helper function to advance time by @p time_step second.
   */
  void evolve_one_time_step(double const time_step,
                            SuperCapacitorState supercapacitor_state,
                            bool rebuild);

  /**
   * Output on the screen the condition number of the system of equations being
   * solved.
   */
  void output_condition_number(double condition_number);

  /**
   * Output on the screen the eigenvalues of the system of equations being
   * solved.
   */
  void output_eigenvalues(std::vector<double> const &eigenvalues);

  /**
   * Helper function for the constructor.
   */
  void setup();

  /**
   * Maximum number of iterations of the Krylov solver in
   * evolve_one_time_step().
   */
  unsigned int _max_iter;
  /**
   * Verbosity level of the Krylov solver in evolve_one_time_step().
   */
  unsigned int _verbose_lvl;
  /**
   * Absolute tolerance of the Krylov solver in evolve_one_time_step(). The
   * tolerance used by the Krylov solver is the maximum of the relative and the
   * absolute tolerance.
   */
  double _abs_tolerance;
  /**
   * Relative tolerance of the Krylov solver in evolve_one_time_step(), i.e. the
   * tolerance is @p rel_tolerance \f$ \times ||b||_{2}\f$. The tolerance used
   * by the Krylov solver is the maximum of the relative and the absolute
   * tolerance.
   */
  double _rel_tolerance;
  /**
   * Area of the cathode.
   */
  double _surface_area;

  std::shared_ptr<Geometry<dim>> _geometry;
  std::shared_ptr<dealii::FESystem<dim>> _fe;
  std::shared_ptr<dealii::DoFHandler<dim>> _dof_handler;
  std::shared_ptr<dealii::Trilinos::MPI::BlockVector> _solution;

  std::shared_ptr<ElectrochemicalPhysicsParameters<dim>>
      _electrochemical_physics_params;
  std::shared_ptr<ElectrochemicalPhysics<dim>> _electrochemical_physics;
  std::shared_ptr<SuperCapacitorPostprocessorParameters<dim>>
      _post_processor_params;
  std::shared_ptr<SuperCapacitorPostprocessor<dim>> _post_processor;
  boost::property_tree::ptree const _ptree;
  Timer _setup_timer;
  Timer _solver_timer;

  template <int dimension>
  friend class SuperCapacitorInspector;
};
}

#endif
