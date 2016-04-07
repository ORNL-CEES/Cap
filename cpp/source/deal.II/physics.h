#ifndef CAP_PHYSICS_H
#define CAP_PHYSICS_H

#include <cap/mp_values.h>
#include <cap/types.h>
#include <boost/property_tree/ptree.hpp>
#include <deal.II/base/index_set.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

namespace cap
{
/**
 * This class encapsulate all the parameters needed to build a Physics object.
 */
template <int dim>
class PhysicsParameters
{
public:
  PhysicsParameters(boost::property_tree::ptree const &d) : database(d) {}

  virtual ~PhysicsParameters() = default;

  std::shared_ptr<Geometry<dim> const> geometry;
  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;
  std::shared_ptr<MPValues<dim> const> mp_values;
  boost::property_tree::ptree const database;
};

/**
 * This is the base class for the Physics. The derived Physics need to build
 * the matrix of the system and the right-hand side. They are also responsible
 * of the boundary conditions.
 */
template <int dim>
class Physics
{
public:
  Physics(std::shared_ptr<PhysicsParameters<dim> const> parameters,
          boost::mpi::communicator mpi_communicator);

  virtual ~Physics() = default;

  inline boost::mpi::communicator get_mpi_communicator() const
  {
    return mpi_communicator;
  }

  inline dealii::Trilinos::SparseMatrix const &get_system_matrix() const
  {
    return system_matrix;
  }

  inline dealii::Trilinos::SparseMatrix const &get_mass_matrix() const
  {
    return mass_matrix;
  }

  inline dealii::ConstraintMatrix const &get_constraint_matrix() const
  {
    return constraint_matrix;
  }

  /**
   * Return the right-hand side of system of equations.
   */
  inline dealii::Trilinos::MPI::Vector const &get_system_rhs() const
  {
    return this->system_rhs;
  }

protected:
  boost::mpi::communicator mpi_communicator;
  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;
  dealii::IndexSet locally_owned_dofs;
  dealii::IndexSet locally_relevant_dofs;
  dealii::ConstraintMatrix constraint_matrix;
  dealii::Trilinos::SparsityPattern sparsity_pattern;
  dealii::Trilinos::SparseMatrix system_matrix;
  dealii::Trilinos::SparseMatrix mass_matrix;
  dealii::Trilinos::MPI::Vector system_rhs;
  std::shared_ptr<MPValues<dim> const> mp_values;
};
}

#endif
