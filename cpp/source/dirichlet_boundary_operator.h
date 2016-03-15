/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <boost/property_tree/ptree.hpp>
#include <map>
#include <memory>

namespace cap
{

class DirichletBoundaryOperator
{
public:
  DirichletBoundaryOperator(
      std::map<dealii::types::global_dof_index, double> const &bv,
      std::shared_ptr<boost::property_tree::ptree> database)
      : boundary_values(bv)
  {
    symmetric_correction = database->get<bool>("symmetric_correction");
    zero_diagonal        = database->get<bool>("zero_diagonal");
  }
  void apply_rhs_correction(dealii::Vector<double> &right_hand_side)
  {
    if (!zero_diagonal)
      for (auto const &x : rhs_set)
        right_hand_side[x.first] = x.second;
    else if (!rhs_set.empty())
      throw dealii::ExcMessage(
          "rhs_set is not empty but zero_diagonal is true");

    if (symmetric_correction)
      for (auto const &x : rhs_add)
        right_hand_side[x.first] += x.second;
    else if (!rhs_add.empty())
      throw dealii::ExcMessage(
          "rhs_add is not empty but symmetric_correction is false");
  }

  void apply_matrix_correction(dealii::SparseMatrix<double> &matrix)
  {
    Assert(matrix.n() == matrix.m(),
           dealii::ExcDimensionMismatch(matrix.n(), matrix.m()));

    if (boundary_values.empty())
      return;

    dealii::types::global_dof_index const n_dofs = matrix.m();

    for (auto const &dof : boundary_values)
    {
      dealii::types::global_dof_index const dof_index = dof.first;
      double const dof_value = dof.second;
      if (!(dof_index < n_dofs))
        throw std::runtime_error("invalid dof index (probably should make this "
                                 "debug only but it helps getting rid of silly "
                                 "warning non used variable)");
      // for each boundary dof:

      // set entries of this line to zero except for the diagonal
      // entry
      for (dealii::SparseMatrix<double>::iterator p = matrix.begin(dof_index);
           p != matrix.end(dof_index); ++p)
        if ((zero_diagonal ? true : (p->column() != dof_index)))
          p->value() = 0.0;

      // set right hand side to
      // wanted value: if main diagonal
      // entry nonzero, don't touch it
      // and scale rhs accordingly. If
      // zero, take the first main
      // diagonal entry we can find, or
      // one if no nonzero main diagonal
      // element exists. Normally, however,
      // the main diagonal entry should
      // not be zero.
      //
      // store the new rhs entry to make
      // the gauss step more efficient
      if (!zero_diagonal)
      {
        double const diagonal_entry = matrix.diag_element(dof_index);
        if (matrix.diag_element(dof_index) != 0.0)
        {
          this->rhs_set[dof_index] = dof_value * diagonal_entry;
        }
        else
        {
          throw dealii::ExcMessage("diagonal entry is zero... "
                                   "not implemented yet");
        }
      }

      // if the user wants to have
      // the symmetry of the matrix
      // preserved, and if the
      // sparsity pattern is
      // symmetric, then do a Gauss
      // elimination step with the
      // present row
      if (this->symmetric_correction)
      {
        // we have to loop over all rows of the matrix which have
        // a nonzero entry in the column which we work in
        // presently. if the sparsity pattern is symmetric, then
        // we can get the positions of these rows cheaply by
        // looking at the nonzero column numbers of the present
        // row. we need not look at the first entry of each row,
        // since that is the diagonal element and thus the present
        // row
        for (dealii::SparseMatrix<double>::iterator q =
                 matrix.begin(dof_index) + 1;
             q != matrix.end(dof_index); ++q)
        {
          const dealii::types::global_dof_index row = q->column();

          // find the position of
          // element
          // (row,dof_double number)
          const dealii::SparseMatrix<double>::iterator p =
              dealii::Utilities::lower_bound(
                  matrix.begin(row) + 1, matrix.end(row), dof_index,
                  [](const dealii::SparseMatrix<double>::iterator::value_type p,
                     const unsigned int column)
                  {
                    return (p.column() < column);
                  });

          // check whether this line has an entry in the
          // regarding column (check for ==dof_number and !=
          // next_row, since if row==dof_number-1, *p is a
          // past-the-end pointer but points to dof_number
          // anyway...)
          //
          // there should be such an entry! we know this because
          // we have assumed that the sparsity pattern is
          // symmetric and we only walk over those rows for
          // which the current row has a column entry
          Assert((p != matrix.end(row)) && (p->column() == dof_index),
                 dealii::ExcMessage(
                     "This function is trying to access an element of the "
                     "matrix that doesn't seem to exist. Are you using a "
                     "nonsymmetric sparsity pattern? If so, you are not "
                     "allowed to set the eliminate_column argument of this "
                     "function, see the documentation."));

          // correct right hand side
          this->rhs_add[row] -= p->value() * dof_value;

          // set matrix entry to zero
          p->value() = 0.0;
        }
      }
    }
  }

private:
  std::map<dealii::types::global_dof_index, double> boundary_values;
  // preserve symmetry
  bool symmetric_correction;
  // zero out the entire row including the diagonal
  bool zero_diagonal;
  std::map<dealii::types::global_dof_index, double> rhs_set;
  std::map<dealii::types::global_dof_index, double> rhs_add;
};

} // end namespace cap
