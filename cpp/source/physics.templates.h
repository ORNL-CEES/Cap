#ifndef CAP_PHYSICS_TEMPLATES_H
#define CAP_PHYSICS_TEMPLATES_H

#include <cap/physics.h>

namespace cap
{
template <int dim>
Physics<dim>::Physics(
    std::shared_ptr<PhysicsParameters<dim> const> parameters)
    : dof_handler(parameters->dof_handler),
      mp_values(parameters->mp_values)
{
  this->system_matrix.reinit(*(parameters->sparsity_pattern));
  this->load_vector.reinit(*(parameters->some_vector));
}
}

#endif
