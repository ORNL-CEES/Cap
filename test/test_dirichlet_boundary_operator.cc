#define BOOST_TEST_MODULE DiricheltBoundaryOperator
#define BOOST_TEST_MAIN
#include <cap/dirichlet_boundary_operator.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <boost/test/unit_test.hpp>

void initialize_system(
    dealii::SparsityPattern      & sparsity_pattern,
    dealii::SparseMatrix<double> & system_matrix, 
    dealii::Vector      <double> & system_rhs   ,
    std::map<dealii::types::global_dof_index,double> & boundary_values
    )
{
    // make grid
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube (triangulation, -1, 1);
    triangulation.refine_global (5);
    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl;

    // setup system
    dealii::FE_Q      <2> fe(1);
    dealii::DoFHandler<2> dof_handler (triangulation);

    dof_handler.distribute_dofs (fe);
    std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;
    dealii::DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
    dealii::DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit (sparsity_pattern);
    system_rhs   .reinit (dof_handler.n_dofs());

    // assemble system
    dealii::QGauss  <2> quadrature_formula(2);
    dealii::FEValues<2> fe_values (fe, quadrature_formula,
                           dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();
    dealii::FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    dealii::Vector    <double>   cell_rhs    (dofs_per_cell);
    std::vector<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);
    dealii::DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;
        cell_rhs = 0;
        for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
          {
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                     fe_values.shape_grad (j, q_index) *
                                     fe_values.JxW (q_index));
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              cell_rhs(i) += (fe_values.shape_value (i, q_index) *
                              1 *
                              fe_values.JxW (q_index));
          }
        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }
    dealii::VectorTools::interpolate_boundary_values (dof_handler,
                                                      0,
                                                      dealii::ConstantFunction<2>(33.0),
//                                                      dealii::ZeroFunction<2>(),
                                                      boundary_values);
}

BOOST_AUTO_TEST_CASE( test_dirichlet_boundary_operator )
{
    dealii::SparsityPattern      sparsity_pattern;
    dealii::SparseMatrix<double> system_matrix;
    dealii::Vector      <double> system_rhs   ;
    std::map<dealii::types::global_dof_index,double> boundary_values;

    initialize_system (sparsity_pattern, system_matrix, system_rhs, boundary_values);

    dealii::Vector<double> solution ( system_rhs.size() );

    std::cout<<"rhs = "<<system_rhs.l2_norm()<<"\n";
    std::cout<<"sol = "<<solution  .l2_norm()<<"\n";

    bool const eliminate_columns = false;

    // deal apply boundary values
    dealii::Vector      <double> dealii_apply_boundary_values_system_rhs(system_rhs);
    dealii::Vector      <double> dealii_apply_boundary_values_solution(system_rhs.size());
    dealii::SparseMatrix<double> dealii_apply_boundary_values_system_matrix(sparsity_pattern);
    dealii_apply_boundary_values_system_matrix.copy_from(system_matrix);
    
    dealii::MatrixTools::apply_boundary_values (boundary_values,
                                                dealii_apply_boundary_values_system_matrix,
                                                dealii_apply_boundary_values_solution,
                                                dealii_apply_boundary_values_system_rhs,
                                                eliminate_columns);

    std::cout<<"deal rhs = "<<dealii_apply_boundary_values_system_rhs.l2_norm()<<"\n";
    std::cout<<"deal sol = "<<dealii_apply_boundary_values_solution  .l2_norm()<<"\n";

    dealii_apply_boundary_values_solution = 1.0;
    dealii_apply_boundary_values_system_matrix.vmult(dealii_apply_boundary_values_system_rhs, dealii_apply_boundary_values_solution);
    std::cout<<"deal apply = "<<dealii_apply_boundary_values_system_rhs.l2_norm()<<"\n";

    // cap dirichlet boundary operator
    dealii::Vector      <double> cap_dirichlet_boundary_operator_system_rhs(system_rhs);
    dealii::Vector      <double> cap_dirichlet_boundary_operator_solution(system_rhs.size());
    dealii::SparseMatrix<double> cap_dirichlet_boundary_operator_system_matrix(sparsity_pattern);
    cap_dirichlet_boundary_operator_system_matrix.copy_from(system_matrix);

    std::shared_ptr<boost::property_tree::ptree> boundary_operator_database =
        std::make_shared<boost::property_tree::ptree>();
    boundary_operator_database->put("symmetric_correction", eliminate_columns);
    boundary_operator_database->put("zero_diagonal"       , false            );
    
    cap::DirichletBoundaryOperator boundary_operator (boundary_values, boundary_operator_database);
    boundary_operator.apply_matrix_correction (cap_dirichlet_boundary_operator_system_matrix);
    boundary_operator.apply_rhs_correction (cap_dirichlet_boundary_operator_system_rhs);

    std::cout<<"cap rhs = "<<cap_dirichlet_boundary_operator_system_rhs.l2_norm()<<"\n";
    std::cout<<"cap sol = "<<cap_dirichlet_boundary_operator_solution  .l2_norm()<<"\n";

    cap_dirichlet_boundary_operator_solution = 1.0;
    cap_dirichlet_boundary_operator_system_matrix.vmult(cap_dirichlet_boundary_operator_system_rhs, cap_dirichlet_boundary_operator_solution);
    std::cout<<"cap apply = "<<cap_dirichlet_boundary_operator_system_rhs.l2_norm()<<"\n";

    double const tolerance = 1.0e-14;
    cap_dirichlet_boundary_operator_system_matrix.add(-1.0, dealii_apply_boundary_values_system_matrix);
    cap_dirichlet_boundary_operator_system_rhs   .add(-1.0, dealii_apply_boundary_values_system_rhs   );
    BOOST_CHECK_SMALL(cap_dirichlet_boundary_operator_system_matrix.l1_norm(), tolerance);
    BOOST_CHECK_SMALL(cap_dirichlet_boundary_operator_system_rhs   .l1_norm(), tolerance);
}
