#ifndef CAP_OPERATOR_H
#define CAP_OPERATOR_H

#include <cap/mp_values.h>
#include <cap/boundary_values.h>

//#include <deal.II/base/types.h>
//#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/lac/constraint_matrix.h>
//#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
//#include <deal.II/lac/vector.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace cap
{

//////////////////////// OPERATOR PARAMETERS ////////////////////////////
template <int dim>
class OperatorParameters
{
public:
  OperatorParameters(std::shared_ptr<boost::property_tree::ptree const> d)
      : database(d)
  {
  }

  virtual ~OperatorParameters() = default;

  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;

  std::shared_ptr<dealii::ConstraintMatrix const> constraint_matrix;

  dealii::SparsityPattern const *sparsity_pattern;

  dealii::Vector<double> const *some_vector;

  std::shared_ptr<MPValues<dim> const> mp_values;

  std::shared_ptr<BoundaryValues<dim> const> boundary_values;

  std::shared_ptr<boost::property_tree::ptree const> database;
};

//////////////////////// OPERATOR ////////////////////////////
/**
 * This is the base class for the operators. The derived operators need to
 * build the mass matrix, the stiffness matrix, and the load vector
 * (right-hand-side vector). These operators are used by EnergyStorageDevice
 * object.
 */
template <int dim>
class Operator
{
public:
  Operator(std::shared_ptr<OperatorParameters<dim> const> parameters);

  virtual ~Operator() = default;

  virtual void reset(std::shared_ptr<OperatorParameters<dim> const>) {}

  inline dealii::SparseMatrix<double> const &get_mass_matrix() const
  {
    return this->mass_matrix;
  }

  inline dealii::SparseMatrix<double> const &get_stiffness_matrix() const
  {
    return this->stiffness_matrix;
  }

  inline dealii::Vector<double> const &get_load_vector() const
  {
    return this->load_vector;
  }

  inline std::map<dealii::types::global_dof_index, double> const &
  get_boundary_values() const
  {
    return this->boundary_values;
  }

  inline std::vector<dealii::types::global_dof_index> const &
  get_null_space() const
  {
    return this->null_space_dof_indices;
  }

  void set_null_space(unsigned int const, dealii::types::material_id const);

protected:
  std::shared_ptr<dealii::DoFHandler<dim> const> dof_handler;
  dealii::types::global_dof_index dof_shift;
  std::shared_ptr<dealii::ConstraintMatrix const> constraint_matrix;

  std::shared_ptr<MPValues<dim> const> mp_values;
  std::shared_ptr<BoundaryValues<dim> const> b_values;

  dealii::SparseMatrix<double> mass_matrix;
  dealii::SparseMatrix<double> stiffness_matrix;
  dealii::Vector<double> load_vector;
  std::map<dealii::types::global_dof_index, double> boundary_values;
  std::vector<dealii::types::global_dof_index> null_space_dof_indices;
};

} // end namespace cap

#endif // CAP_OPERATOR_H
