#ifndef CAP_DEAL_II_ELECTROCHEMICAL_PHYSICS_H
#define CAP_DEAL_II_ELECTROCHEMICAL_PHYSICS_H

#include <cap/physics.h>

namespace cap
{
  /**
   * This class builds the system of equations that describes an electrochemical
   * physics. The system is built when the constructor or the reinit() function
   * is called.
   */
  template <int dim>
  class ElectrochemicalPhysics : public Physics<dim>
  {
    public:
      ElectrochemicalPhysics(std::shared_ptr<PhysicsParameters<dim> const> parameters);

    private:
      void assemble_system();

      unsigned int solid_potential_component;
      unsigned int liquid_potential_component;
      unsigned int temperature_component;
      double alpha;
      dealii::types::boundary_id anode_boundary_id;
      dealii::types::boundary_id cathode_boundary_id;
  };
}

#endif
