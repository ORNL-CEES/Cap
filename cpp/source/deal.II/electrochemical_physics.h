/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_DEAL_II_ELECTROCHEMICAL_PHYSICS_H
#define CAP_DEAL_II_ELECTROCHEMICAL_PHYSICS_H

#include <cap/physics.h>
#include <cap/timer.h>

namespace cap
{
enum SuperCapacitorState
{
  Uninitialized,
  ConstantCurrent,
  ConstantVoltage,
  ConstantLoad
};

/**
 * This class encapsulates the parameters used in ElectrochemicalPhysics.
 */
template <int dim>
class ElectrochemicalPhysicsParameters : public PhysicsParameters<dim>
{
public:
  ElectrochemicalPhysicsParameters(boost::property_tree::ptree const &d)
      : PhysicsParameters<dim>(d), supercapacitor_state(Uninitialized),
        constant_current_density(0.), constant_voltage(0.),
        constant_load_density(0.), time_step(0.)
  {
  }

  SuperCapacitorState supercapacitor_state;
  double constant_current_density;
  double constant_voltage;
  double constant_load_density;
  double time_step;
};

/**
 * This class builds the system of equations that describes an electrochemical
 * physics. The system is built when the constructor or the reinit() function
 * is called.
 */
template <int dim>
class ElectrochemicalPhysics : public Physics<dim>
{
public:
  ElectrochemicalPhysics(
      std::shared_ptr<PhysicsParameters<dim> const> parameters,
      boost::mpi::communicator mpi_communicator);

  ~ElectrochemicalPhysics();

private:
  void assemble_system(std::shared_ptr<PhysicsParameters<dim> const> parameters,
                       bool const inhomogeneous_bc);

  unsigned int solid_potential_component;
  unsigned int liquid_potential_component;
  Timer _assembly_timer;
  Timer _setup_timer;
};
}

#endif
