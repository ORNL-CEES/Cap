#include <cap/dof_extractor.templates.h>

namespace cap
{

DoFExtractor::DoFExtractor(dealii::ComponentMask const &row_mask,
                           dealii::ComponentMask const &column_mask,
                           unsigned int const dofs_per_cell)
{
  this->reset(row_mask, column_mask, dofs_per_cell);
}

void DoFExtractor::reset(dealii::ComponentMask const &row_mask,
                         dealii::ComponentMask const &column_mask,
                         unsigned int const dofs_per_cell)
{
  unsigned int const n_components = row_mask.size();
  Assert(column_mask.size() == n_components,
         dealii::StandardExceptions::ExcDimensionMismatch(row_mask.size(),
                                                          column_mask.size()));
  this->row_indices.clear();
  this->row_indices.reserve(dofs_per_cell);
  this->column_indices.clear();
  this->column_indices.reserve(dofs_per_cell);
  for (unsigned int dof = 0; dof < dofs_per_cell; ++dof)
  {
    if (row_mask[dof % n_components])
    {
      this->row_indices.push_back(dof);
    } // end if
    if (column_mask[dof % n_components])
    {
      this->column_indices.push_back(dof);
    } // end if
  }   // end for dof
}

std::vector<dealii::types::global_dof_index> DoFExtractor::extract_row_indices(
    std::vector<dealii::types::global_dof_index> const &dofs_in) const
{
  unsigned int const n_dofs = this->row_indices.size();
  std::vector<dealii::types::global_dof_index> dofs_out(n_dofs);
  for (unsigned int dof = 0; dof < n_dofs; ++dof)
  {
    dofs_out[dof] = dofs_in[this->row_indices[dof]];
  } // end for dof
  return dofs_out;
}

std::vector<dealii::types::global_dof_index>
DoFExtractor::extract_column_indices(
    std::vector<dealii::types::global_dof_index> const &dofs_in) const
{
  unsigned int const n_dofs = this->column_indices.size();
  std::vector<dealii::types::global_dof_index> dofs_out(n_dofs);
  for (unsigned int dof = 0; dof < n_dofs; ++dof)
  {
    dofs_out[dof] = dofs_in[this->column_indices[dof]];
  } // end for dof
  return dofs_out;
}

template dealii::FullMatrix<double>
DoFExtractor::extract_matrix(dealii::FullMatrix<double> const &) const;
template dealii::FullMatrix<float>
DoFExtractor::extract_matrix(dealii::FullMatrix<float> const &) const;
template dealii::Vector<double>
DoFExtractor::extract_vector(dealii::Vector<double> const &) const;
template dealii::Vector<float>
DoFExtractor::extract_vector(dealii::Vector<float> const &) const;

} // end namespace cap
