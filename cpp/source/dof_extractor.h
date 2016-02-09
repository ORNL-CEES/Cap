/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_DOF_EXTRACTOR_H
#define CAP_DOF_EXTRACTOR_H

#include <deal.II/base/types.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

namespace cap
{

class DoFExtractor
{
public:
  DoFExtractor(dealii::ComponentMask const &row_mask,
               dealii::ComponentMask const &column_mask,
               unsigned int const dofs_per_cell);
  void reset(dealii::ComponentMask const &row_mask,
             dealii::ComponentMask const &column_mask,
             unsigned int const dofs_per_cell);
  std::vector<dealii::types::global_dof_index> extract_row_indices(
      std::vector<dealii::types::global_dof_index> const &) const;
  std::vector<dealii::types::global_dof_index> extract_column_indices(
      std::vector<dealii::types::global_dof_index> const &) const;
  template <typename T>
  dealii::FullMatrix<T> extract_matrix(dealii::FullMatrix<T> const &) const;
  template <typename T>
  dealii::Vector<T> extract_vector(dealii::Vector<T> const &) const;

private:
  std::vector<unsigned int> row_indices;
  std::vector<unsigned int> column_indices;
};

} // end namespace cap

#endif // CAP_DOF_EXTRACTOR_H
