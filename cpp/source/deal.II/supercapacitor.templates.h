/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_DEAL_II_SUPERCAPACITOR_TEMPLATES_H
#define CAP_DEAL_II_SUPERCAPACITOR_TEMPLATES_H

#include <cap/supercapacitor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>
#include <fstream>
#include <typeinfo>

namespace cap
{
template <int dim>
void SuperCapacitorInspector<dim>::inspect(EnergyStorageDevice *device)
{
  static int i = 0;
  SuperCapacitor<dim> *supercapacitor =
      dynamic_cast<SuperCapacitor<dim> *>(device);
  // dynamic_cast does not throw an exception when using pointer. It just sets
  // the pointer to nullptr, so we throw the bad_cast exception ourselves.
  if (supercapacitor == nullptr)
    throw std::bad_cast();

  std::vector<std::string> keys =
      supercapacitor->post_processor->get_vector_keys();
  std::shared_ptr<dealii::distributed::Triangulation<dim> const> triangulation =
      supercapacitor->_geometry->get_triangulation();
  dealii::DataOut<dim> data_out;
  data_out.attach_triangulation(*triangulation);
  // Output the subdomain id
  dealii::Vector<float> subdomain(triangulation->n_active_cells());
  dealii::types::subdomain_id const local_subdomain_id =
      triangulation->locally_owned_subdomain();
  for (auto &subdom : subdomain)
    subdom = local_subdomain_id;
  data_out.add_data_vector(subdomain, "subdomain");
  // Output the required quantities
  if (!keys.empty())
  {
    BOOST_FOREACH (std::string const &key, keys)
      data_out.add_data_vector(supercapacitor->post_processor->get(key), key);
  }
  data_out.build_patches();
  std::string const filename =
      "solution-" + dealii::Utilities::int_to_string(i, 4) + "." +
      dealii::Utilities::int_to_string(local_subdomain_id, 4) + ".vtu";
  std::ofstream fout(filename.c_str());
  data_out.write_vtu(fout);
  fout.close();

  // Create the master record
  if (supercapacitor->_communicator.rank() == 0)
  {
    std::vector<std::string> filenames;
    for (int j = 0; j < supercapacitor->_communicator.size(); ++j)
      filenames.push_back("solution-" + dealii::Utilities::int_to_string(i, 4) +
                          "." + dealii::Utilities::int_to_string(j, 4) +
                          ".vtu");
    std::ofstream master_output(
        ("solution-" + dealii::Utilities::int_to_string(i, 4) + ".pvtu")
            .c_str());
    data_out.write_pvtu_record(master_output, filenames);
  }
  ++i;
}

template <int dim>
SuperCapacitor<dim>::SuperCapacitor(boost::property_tree::ptree const &ptree,
                                    boost::mpi::communicator const &comm)
    : EnergyStorageDevice(comm), max_iter(0), verbose_lvl(0), abs_tolerance(0.),
      rel_tolerance(0.), surface_area(0.), _geometry(nullptr), _fe(nullptr),
      dof_handler(nullptr), solution(nullptr),
      electrochemical_physics_params(nullptr), electrochemical_physics(nullptr),
      post_processor_params(nullptr), post_processor(nullptr), _ptree(ptree),
      _setup_timer(comm, "SuperCapacitor setup"),
      _solver_timer(comm, "SuperCapacitor solver")
{
  _setup_timer.start();

  verbose_lvl = _ptree.get("verbosity", 0);

  // get data tolerance and maximum number of iterations for the CG solver
  boost::property_tree::ptree const &solver_database =
      _ptree.get_child("solver");
  max_iter = solver_database.get("max_iter", 1000);
  rel_tolerance = solver_database.get("rel_tolerance", 1e-12);
  abs_tolerance = solver_database.get("abs_tolerance", 1e-12);
  // set the number of threads used by deal.II
  unsigned int n_threads = solver_database.get("n_threads", 1);
  // if 0, let TBB uses all the available threads. This can also be used if one
  // wants to use DEAL_II_NUM_THREADS
  if (n_threads == 0)
    dealii::MultithreadInfo::set_thread_limit();
  else
    dealii::MultithreadInfo::set_thread_limit(n_threads);

  // build triangulation
  std::shared_ptr<boost::property_tree::ptree> geometry_database =
      std::make_shared<boost::property_tree::ptree>(
          _ptree.get_child("geometry"));
  _geometry = std::make_shared<cap::Geometry<dim>>(geometry_database,
                                                   this->_communicator);
  std::string mesh_type = geometry_database->get<std::string>("type");
  if (mesh_type.compare("restart") != 0)
    setup();
}

template <int dim>
SuperCapacitor<dim>::~SuperCapacitor()
{
  if (verbose_lvl > 0)
  {
    _setup_timer.print();
    _solver_timer.print();
  }
}

template <int dim>
void SuperCapacitor<dim>::inspect(EnergyStorageDeviceInspector *inspector)
{
  inspector->inspect(this);
}

template <int dim>
void SuperCapacitor<dim>::get_voltage(double &voltage) const
{
  post_processor->get("voltage", voltage);
}

template <int dim>
void SuperCapacitor<dim>::get_current(double &current) const
{
  post_processor->get("current", current);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_constant_current(
    double const time_step, double const current)
{
  BOOST_ASSERT_MSG(surface_area > 0.,
                   "The surface area should be greater than zero.");
  double const constant_current_density = current / surface_area;
  bool const rebuild =
      (electrochemical_physics_params->constant_current_density ==
       constant_current_density)
          ? false
          : true;
  electrochemical_physics_params->constant_current_density =
      constant_current_density;
  evolve_one_time_step(time_step, ConstantCurrent, rebuild);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_constant_voltage(
    double const time_step, double const voltage)
{
  bool const rebuild =
      (electrochemical_physics_params->constant_voltage == voltage) ? false
                                                                    : true;
  electrochemical_physics_params->constant_voltage = voltage;
  evolve_one_time_step(time_step, ConstantVoltage, rebuild);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_constant_power(
    double const time_step, double const power)
{
  BOOST_ASSERT_MSG(surface_area > 0.,
                   "The surface area should be greater than zero.");
  dealii::Trilinos::MPI::Vector old_solution(solution->block(0));
  // The tolerance and the maximum number of iterations are for the picard
  // iterations done below. This is not related to the Krylov solver in
  // evolve_one_time_step.
  int const max_iterations = 10;
  double const percent_tolerance = 1.0e-2;
  double current(0.0);
  double voltage(0.0);
  get_voltage(voltage);
  for (int k = 0; k < max_iterations; ++k)
  {
    current = power / voltage;
    double const constant_current_density = current / surface_area;
    bool const rebuild =
        (electrochemical_physics_params->constant_current_density ==
         constant_current_density)
            ? false
            : true;
    electrochemical_physics_params->constant_current_density =
        constant_current_density;
    evolve_one_time_step(time_step, ConstantCurrent, rebuild);
    get_voltage(voltage);
    if (std::abs(power - voltage * current) / std::abs(power) <
        percent_tolerance)
      return;
    solution->block(0) = old_solution;
  }
  throw std::runtime_error("fixed point iteration did not converge in " +
                           std::to_string(max_iterations) + " iterations");
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_constant_load(
    double const time_step, double const load)
{
  electrochemical_physics_params->constant_load_density = load * surface_area;
  // BC not implemented yet
  throw std::runtime_error("This function is not implemented.");
  evolve_one_time_step(time_step, ConstantLoad, true);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_linear_current(
    double const time_step, double const current)
{
  // TODO: this is a temporary solution
  evolve_one_time_step_constant_current(time_step, current);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_linear_voltage(
    double const time_step, double const voltage)
{
  // TODO: this is a temporary solution
  evolve_one_time_step_constant_voltage(time_step, voltage);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_linear_power(
    double const time_step, double const power)
{
  // TODO: this is a temporary solution
  evolve_one_time_step_constant_power(time_step, power);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step_linear_load(
    double const time_step, double const load)
{
  // TODO: this is a temporary solution
  evolve_one_time_step_constant_load(time_step, load);
}

template <int dim>
void SuperCapacitor<dim>::evolve_one_time_step(
    double const time_step, SuperCapacitorState supercapacitor_state,
    bool rebuild)
{
  // The first time evolve_one_time_step is called, the solution and the
  // post-processor need to be iniatialized.
  if (electrochemical_physics_params->supercapacitor_state == Uninitialized)
  {
    electrochemical_physics_params->time_step = time_step;
    electrochemical_physics_params->supercapacitor_state = supercapacitor_state;
    electrochemical_physics.reset(new ElectrochemicalPhysics<dim>(
        electrochemical_physics_params, this->_communicator));
  }
  // Rebuild the system if necessary
  else if ((rebuild == true) ||
           (std::abs(time_step / electrochemical_physics_params->time_step -
                     1.0) > 1e-14) ||
           (supercapacitor_state !=
            electrochemical_physics_params->supercapacitor_state))
  {
    electrochemical_physics_params->time_step = time_step;
    electrochemical_physics_params->supercapacitor_state = supercapacitor_state;
    electrochemical_physics.reset(new ElectrochemicalPhysics<dim>(
        electrochemical_physics_params, this->_communicator));
  }

  // Get the system from the ElectrochemicalPhysiscs object.
  dealii::Trilinos::SparseMatrix const &system_matrix =
      electrochemical_physics->get_system_matrix();
  dealii::Trilinos::SparseMatrix const &mass_matrix =
      electrochemical_physics->get_mass_matrix();
  dealii::ConstraintMatrix const &constraint_matrix =
      electrochemical_physics->get_constraint_matrix();
  dealii::Trilinos::MPI::Vector const &system_rhs =
      electrochemical_physics->get_system_rhs();
  dealii::Trilinos::MPI::Vector time_dep_rhs = system_rhs;
  mass_matrix.vmult_add(time_dep_rhs, solution->block(0));

  // Solve the system
  _solver_timer.start();
  double tolerance =
      std::max(abs_tolerance, rel_tolerance * system_rhs.l2_norm());
  dealii::SolverControl solver_control(max_iter, tolerance);
  dealii::SolverCG<dealii::Trilinos::MPI::Vector> solver(solver_control);
  // Compute the condition number at the end of the CG iterations.
  if (verbose_lvl > 1)
    solver.connect_condition_number_slot(
        std::bind(&SuperCapacitor<dim>::output_condition_number, this,
                  std::placeholders::_1),
        false);
  // Compute all the eigenvalues at the end of the CG iterations.
  if (verbose_lvl > 2)
    solver.connect_eigenvalues_slot(
        std::bind(&SuperCapacitor<dim>::output_eigenvalues, this,
                  std::placeholders::_1),
        false);
  dealii::Trilinos::PreconditionAMG preconditioner;
  // Temporary preconditioner. Need to find what parameters work best.
  preconditioner.initialize(system_matrix);
  constraint_matrix.distribute(solution->block(0));
  solver.solve(system_matrix, solution->block(0), time_dep_rhs, preconditioner);
  constraint_matrix.distribute(solution->block(0));
  if ((verbose_lvl > 0) && (_communicator.rank() == 0))
  {
    std::cout << "Initial value: " << solver_control.initial_value()
              << std::endl;
    std::cout << "Last value: " << solver_control.last_value() << std::endl;
    std::cout << "Number of iterations: " << solver_control.last_step()
              << std::endl
              << std::endl;
  }
  _solver_timer.stop();

  // Update the data in post-processor
  post_processor->reset(post_processor_params);
}

template <int dim>
void SuperCapacitor<dim>::output_condition_number(double condition_number)
{
  if (_communicator.rank() == 0)
    std::cout << "Condition number: " << condition_number << std::endl;
}

template <int dim>
void SuperCapacitor<dim>::output_eigenvalues(
    std::vector<double> const &eigenvalues)
{
  if (_communicator.rank() == 0)
  {
    std::cout << "Eigenvalues:";
    for (auto eig : eigenvalues)
      std::cout << " " << eig;
    std::cout << std::endl;
  }
}

template <int dim>
std::shared_ptr<Geometry<dim>> SuperCapacitor<dim>::get_geometry() const
{
  return _geometry;
}

template <int dim>
std::shared_ptr<PostprocessorParameters<dim>>
SuperCapacitor<dim>::get_post_processor_parameters() const
{
  return post_processor_params;
}

template <int dim>
std::shared_ptr<Postprocessor<dim>>
SuperCapacitor<dim>::get_post_processor() const
{
  return post_processor;
}

template <int dim>
boost::property_tree::ptree const *
SuperCapacitor<dim>::get_property_tree() const
{
  return &_ptree;
}

template <int dim>
void SuperCapacitor<dim>::save(const std::string &filename) const
{
  // This only stores the refinement information and the solution. The
  // triangulation needs to exist and to be of cells of only the coarsest level,
  // i.e., the coarsest mesh as possible. The vectors that are serialized need
  // to be ghosted.
  unsigned int const n_blocks = solution->n_blocks();
  dealii::IndexSet locally_owned_dofs = dof_handler->locally_owned_dofs();
  dealii::IndexSet locally_relevant_dofs;
  dealii::DoFTools::extract_locally_relevant_dofs(*dof_handler,
                                                  locally_relevant_dofs);
  std::vector<dealii::IndexSet> locally_owned_index_sets(n_blocks,
                                                         locally_owned_dofs);
  std::vector<dealii::IndexSet> locally_relevant_index_sets(
      n_blocks, locally_relevant_dofs);
  dealii::Trilinos::MPI::BlockVector ghosted_solution(
      locally_owned_index_sets, locally_relevant_index_sets,
      this->_communicator);
  ghosted_solution = *solution;

  dealii::distributed::SolutionTransfer<dim, dealii::Trilinos::MPI::BlockVector>
      solution_transfer(*dof_handler);
  solution_transfer.prepare_serialization(ghosted_solution);
  _geometry->get_triangulation()->save(filename.c_str());
}

template <int dim>
void SuperCapacitor<dim>::load(const std::string &filename)
{
  namespace boost_io = boost::iostreams;

  // Check that the files exist
  std::string const coarse_mesh_filename =
      _ptree.get<std::string>("geometry.coarse_mesh_filename");
  if (boost::filesystem::exists(coarse_mesh_filename) == false)
    throw std::runtime_error("The file " + coarse_mesh_filename +
                             " does not exist.");
  if (boost::filesystem::exists(filename) == false)
    throw std::runtime_error("The file " + filename + " does not exist.");

  // Open the file
  std::ifstream is(coarse_mesh_filename, std::ios::binary);
  if (is.good() == false)
    throw std::runtime_error("Error while opening the file: " + filename);

  // Load the coarse mesh. Because the p4est objects are not serialized, we
  // deserialize in a dealii::Triangulation and then use copy_triangulation to
  // copy the dealii::Triangulation into a
  // dealii::parallel::distributed::Triangulation.
  boost_io::filtering_streambuf<boost_io::input> compressed_in;
  compressed_in.push(boost_io::zlib_decompressor());
  compressed_in.push(is);
  boost::archive::binary_iarchive ia(compressed_in);
  dealii::Triangulation<dim> tmp;
  ia >> tmp;
  _geometry->get_triangulation()->clear();
  _geometry->get_triangulation()->copy_triangulation(tmp);

  std::shared_ptr<std::unordered_map<
      std::string, std::set<dealii::types::material_id>>> materials;
  std::shared_ptr<std::unordered_map<
      std::string, std::set<dealii::types::boundary_id>>> boundaries;
  ia >> materials;
  ia >> boundaries;
  _geometry->set_materials(materials);
  _geometry->set_boundaries(boundaries);

  // Load the refinement
  _geometry->get_triangulation()->load(filename.c_str());
  // Do the load balancing. We may want to use different weights after a restart
  // so the weights are not save.
  std::shared_ptr<boost::property_tree::ptree> geometry_database =
      std::make_shared<boost::property_tree::ptree>(
          _ptree.get_child("geometry"));
  _geometry->repartition();

  // Setup the SuperCapacitor object.
  setup();

  // Load the solution.
  dealii::distributed::SolutionTransfer<dim, dealii::Trilinos::MPI::BlockVector>
      solution_transfer(*dof_handler);
  solution_transfer.deserialize(*solution);

  // Reset the post-processor otherwise if we call get_current/get_voltage just
  // after a load the value would stay at zero.
  post_processor->reset(post_processor_params);
}

template <int dim>
void SuperCapacitor<dim>::setup()
{
  std::shared_ptr<dealii::distributed::Triangulation<dim> const> triangulation =
      _geometry->get_triangulation();

  // distribute degrees of freedom
  _fe = std::make_shared<dealii::FESystem<dim>>(dealii::FE_Q<dim>(1), 2);
  dof_handler = std::make_shared<dealii::DoFHandler<dim>>(*triangulation);
  dof_handler->distribute_dofs(*_fe);

  // Renumber the degrees of freedom component-wise.
  dealii::DoFRenumbering::component_wise(*dof_handler);
  unsigned int const n_components =
      dealii::DoFTools::n_components(*dof_handler);
  std::vector<dealii::types::global_dof_index> dofs_per_component(n_components);
  dealii::DoFTools::count_dofs_per_component(*dof_handler, dofs_per_component);

  // read material properties
  std::shared_ptr<boost::property_tree::ptree> material_properties_database =
      std::make_shared<boost::property_tree::ptree>(
          _ptree.get_child("material_properties"));
  MPValuesParameters<dim> params(material_properties_database);
  params.geometry = _geometry;
  std::shared_ptr<MPValues<dim>> mp_values =
      SuperCapacitorMPValuesFactory<dim>::build(params);

  // Initialize the electrochemical physics parameters
  electrochemical_physics_params.reset(
      new ElectrochemicalPhysicsParameters<dim>(_ptree));
  electrochemical_physics_params->geometry = _geometry;
  electrochemical_physics_params->dof_handler = dof_handler;
  electrochemical_physics_params->mp_values =
      std::dynamic_pointer_cast<MPValues<dim> const>(mp_values);

  // Compute the surface area. This is neeeded by several evolve_one_time_step_*
  surface_area = 0.;
  auto const &cathode_boundary_ids = (*_geometry->get_boundaries())["cathode"];
  dealii::QGauss<dim - 1> face_quadrature_rule(_fe->degree + 1);
  unsigned int const n_face_q_points = face_quadrature_rule.size();
  dealii::FEFaceValues<dim> fe_face_values(*_fe, face_quadrature_rule,
                                           dealii::update_JxW_values);
  // TODO this can be simplified when using the next version of deal (current is
  // 8.4)
  for (auto cell : dof_handler->active_cell_iterators())
    if (cell->is_locally_owned() && cell->at_boundary())
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell; ++face)
        if ((cell->face(face)->at_boundary()) &&
            (cathode_boundary_ids.count(cell->face(face)->boundary_id()) > 0))
          for (unsigned int face_q_point = 0; face_q_point < n_face_q_points;
               ++face_q_point)
          {
            fe_face_values.reinit(cell, face);
            surface_area += fe_face_values.JxW(face_q_point);
          }
  // Reduce the value computed on each processor.
  surface_area = dealii::Utilities::MPI::sum(surface_area, this->_communicator);

  // Create the post-processor parameters
  post_processor_params =
      std::make_shared<SuperCapacitorPostprocessorParameters<dim>>(
          std::make_shared<boost::property_tree::ptree>(_ptree), dof_handler);

  // Initialize the size solution
  // Temporary keep using a BlockVector because of PostProcessor.
  std::vector<dealii::IndexSet> index_set(
      1, this->dof_handler->locally_owned_dofs());
  solution.reset(new dealii::Trilinos::MPI::BlockVector(index_set));

  // Initialize the postprocessor
  post_processor_params->solution = solution;
  post_processor_params->mp_values = electrochemical_physics_params->mp_values;
  post_processor = std::make_shared<SuperCapacitorPostprocessor<dim>>(
      post_processor_params, _geometry, this->_communicator);

  post_processor->reset(post_processor_params);

  _setup_timer.stop();
}

} // end namespace cap

#endif
