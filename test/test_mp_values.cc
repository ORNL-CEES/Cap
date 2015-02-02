#define BOOST_TEST_MODULE MaterialPropertyValues
#define BOOST_TEST_DYN_LINK
#include <cap/mp_values.h>
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>

BOOST_AUTO_TEST_CASE( test_mp_values )
{
   std::shared_ptr<boost::property_tree::ptree> empty_database;
   std::shared_ptr<cap::MPValues<2> > mp_values = 
       std::make_shared<cap::MPValues<2> >(empty_database);

   dealii::Triangulation<2> triangulation;
   dealii::GridGenerator::hyper_cube (triangulation);
   dealii::DoFHandler<2> dof_handler(triangulation);
   dealii::DoFHandler<2>::active_cell_iterator cell = 
       dof_handler.begin_active();

   std::vector<double> values;
   std::vector<dealii::Tensor<1, 2> > vectors;
   BOOST_CHECK_THROW( mp_values->get_values("key", cell, values ), std::runtime_error );
   BOOST_CHECK_THROW( mp_values->get_values("key", cell, vectors), std::runtime_error );

   BOOST_CHECK_THROW( 
       std::make_shared<cap::SuperCapacitorMPValues<2> >(empty_database),
       std::runtime_error );

   std::shared_ptr<boost::property_tree::ptree> database(new boost::property_tree::ptree);
   BOOST_CHECK_THROW(
       std::make_shared<cap::SuperCapacitorMPValues<2> >(database),
       boost::property_tree::ptree_bad_path );

}

