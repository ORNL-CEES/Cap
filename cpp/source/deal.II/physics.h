#ifndef CAP_PHYSICS_H
#define CAP_PHYSICS_H

#include <cap/boundary_values.h>
#include <cap/mp_values.h>
#include <boost/property_tree/ptree.hpp>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

namespace cap
{
/**
 * This class encapsulate all the parameters needed to build a Physics object.
 */
template <int dim>
class PhysicsParameters
{
public:
  PhysicsParameters(std::shared_ptr<boost::property_tree::ptree const> d)
      : database(d)
  {
  }

  virtual ~PhysicsParameters() = default;

  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;
  dealii::types::global_dof_index n_dofs;
  std::shared_ptr<MPValues<dim> const> mp_values;
  std::shared_ptr<BoundaryValues<dim> const> boundary_values;
  std::shared_ptr<boost::property_tree::ptree const> database;
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
  Physics(std::shared_ptr<PhysicsParameters<dim> const> parameters);

  virtual ~Physics() = default;

  /**
   * This function is function is not implemented and will throw an
   * exception.
   */
  virtual void reset(std::shared_ptr<PhysicsParameters<dim> const> parameters)
  {
    std::ignore = parameters;
    std::runtime_error("This is function is not implemented");
  }

  inline dealii::SparseMatrix<double> const &get_system_matrix() const
  {
    return system_matrix;
  }

  /**
   * Return the dealii::LinearOperator associated to system_matrix.
   */
  inline dealii::LinearOperator<dealii::Vector<double>, dealii::Vector<double>>
  get_linear_operator() const
  {
    return linear_operator(system_matrix);
  }

  inline dealii::SparseMatrix<double> const &get_mass_matrix() const
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
  inline dealii::Vector<double> const &get_system_rhs() const
  {
    return this->system_rhs;;
  }

protected:
  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;
  dealii::types::global_dof_index n_dofs;
  dealii::ConstraintMatrix constraint_matrix;
  dealii::SparsityPattern sparsity_pattern;
  dealii::SparseMatrix<double> system_matrix;
  dealii::SparseMatrix<double> mass_matrix;
  dealii::Vector<double> system_rhs;
  std::shared_ptr<MPValues<dim> const> mp_values;
  std::shared_ptr<BoundaryValues<dim> const> boundary_values;
};
}

#endif
