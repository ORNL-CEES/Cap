#ifndef CAP_PHYSICS_H
#define CAP_PHYSICS_H

#include <cap/mp_values.h>
#include <boost/property_tree/ptree.hpp>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/sparse_matrix.h>
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
  dealii::SparsityPattern const *sparsity_pattern;
  dealii::Vector<double> const *some_vector;
  std::shared_ptr<MPValues<dim> const> mp_values;
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
    return this->system_matrix;
  }

  /**
   * Return the dealii::LinearOperator associated to system_matrix.
   */
  inline dealii::LinearOperator<dealii::Vector<double>, dealii::Vector<double>>
  get_linear_operator() const
  {
    return linear_operator(this->system_matrix);
  }

  /**
   * Return the right-hand side of system of equations.
   */
  inline dealii::Vector<double> const &get_load_vector() const
  {
    return this->load_vector;
  }

protected:
  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;
  std::shared_ptr<MPValues<dim> const> mp_values;
  dealii::SparseMatrix<double> system_matrix;
  dealii::Vector<double> load_vector;
};
}

#endif
