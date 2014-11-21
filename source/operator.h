#ifndef CACHE_OPERATOR_H
#define CACHE_OPERATOR_H

#include <cache/mp_values.h>
#include <cache/boundary_values.h>

//#include <deal.II/base/types.h>
//#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/lac/constraint_matrix.h>
//#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
//#include <deal.II/lac/vector.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>

//////////////////////// OPERATOR PARAMETERS ////////////////////////////
template <int dim>
class OperatorParameters { 
public:
    OperatorParameters(std::shared_ptr<boost::property_tree::ptree const> d)
        : database(d)
    {  }
    virtual ~OperatorParameters() { }

    dealii::DoFHandler<dim> const * dof_handler;
    dealii::ConstraintMatrix const * constraint_matrix;
    dealii::SparsityPattern const * sparsity_pattern;
    dealii::Vector<double> const * some_vector;

    std::shared_ptr<MPValues<dim> const>       mp_values;
    std::shared_ptr<BoundaryValues<dim> const> boundary_values;

    std::shared_ptr<boost::property_tree::ptree const> database;
};

//////////////////////// OPERATOR ////////////////////////////
template <int dim>
class Operator {
public:
    Operator(OperatorParameters<dim> const & parameters);
    virtual ~Operator() { }
    virtual void reset(OperatorParameters<dim> const & ) { }

    inline dealii::SparseMatrix<double> const & get_mass_matrix() const { return this->mass_matrix; }
    inline dealii::SparseMatrix<double> const & get_stiffness_matrix() const { return this->stiffness_matrix; }
    inline dealii::Vector<double> const & get_load_vector() const { return this->load_vector; }
    inline std::map<dealii::types::global_dof_index, double> const & get_boundary_values() const { return this->boundary_values; }
    inline std::vector<dealii::types::global_dof_index> const & get_null_space() const { return this->null_space_dof_indices; }
    void set_null_space(unsigned int const , dealii::types::material_id const );
    

protected:
    dealii::DoFHandler<dim> const & dof_handler;
    dealii::ConstraintMatrix const & constraint_matrix;

    std::shared_ptr<MPValues<dim> const>       mp_values;
    std::shared_ptr<BoundaryValues<dim> const> b_values;

    dealii::SparseMatrix<double> mass_matrix;
    dealii::SparseMatrix<double> stiffness_matrix;
    dealii::Vector<double> load_vector;
    std::map<dealii::types::global_dof_index, double> boundary_values;
    std::vector<dealii::types::global_dof_index> null_space_dof_indices;
};

#endif // CACHE_OPERATOR_H
