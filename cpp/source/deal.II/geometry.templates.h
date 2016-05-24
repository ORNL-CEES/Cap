/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/geometry.h>
#include <cap/types.h>
#include <cap/utils.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/geometry_info.h>
#include <algorithm>
#include <fstream>
#include <tuple>

namespace cap
{

namespace internal
{
template <int dim>
struct Component
{
public:
  Component(MPI_Comm mpi_communicator)
      : offset(0.), box_dimensions(0), repetitions(0),
        triangulation(mpi_communicator), shift_vector()
  {
  }

  Component(std::vector<dealii::Point<dim>> const &box,
            std::vector<unsigned int> const &repetitions,
            MPI_Comm mpi_communicator)
      : offset(0.), box_dimensions(box), repetitions(repetitions),
        triangulation(mpi_communicator), shift_vector()
  {
  }

  // The current offset of the component.
  double offset;
  std::vector<dealii::Point<dim>> box_dimensions;
  std::vector<unsigned int> repetitions;
  dealii::distributed::Triangulation<dim> triangulation;
  // The last shift that has been applied to the component.
  dealii::Tensor<1, dim> shift_vector;
};

template <int dim>
dealii::Point<dim> transform_coll_a(dealii::Point<dim> const &in,
                                    double const scale_factor,
                                    double const max_value)
{
  double const eps = 1e-15;
  if ((in[dim - 1] < eps) || (in[dim - 1] > (max_value - eps)))
    return in;
  else
  {
    dealii::Point<dim> out = in;
    out[dim - 1] *= scale_factor;
    return out;
  }
}

template <int dim>
dealii::Point<dim> transform_coll_c(dealii::Point<dim> const &in,
                                    double const scale_factor,
                                    double const max_value, double const offset)
{
  double const eps = 1e-15;
  if ((in[dim - 1] < eps) || (in[dim - 1] > (max_value - eps)))
    return in;
  else
  {
    dealii::Point<dim> out = in;
    out[dim - 1] *= scale_factor;
    out[dim - 1] += offset;
    return out;
  }
}

template <int dim>
void read_component_database(boost::property_tree::ptree const &database,
                             Component<dim> &component)
{
  std::vector<unsigned int> repetitions =
      to_vector<unsigned int>(database.get<std::string>("divisions"));
  component.repetitions = repetitions;
  // Add the origin points
  component.box_dimensions.push_back(dealii::Point<dim>());
  std::vector<double> box_dimensions =
      to_vector<double>(database.get<std::string>("dimensions"));
  if (dim == 2)
    component.box_dimensions.push_back(
        dealii::Point<dim>(box_dimensions[0], box_dimensions[1]));
  else
    component.box_dimensions.push_back(dealii::Point<dim>(
        box_dimensions[0], box_dimensions[1], box_dimensions[2]));
}

template <int dim>
void merge_components(
    Component<dim> &component, double const offset,
    std::shared_ptr<dealii::distributed::Triangulation<dim>> _triangulation)
{
  // Move the component by the necessary, i.e. desired position minus the
  // current position
  component.shift_vector[0] = offset - component.offset;
  dealii::GridTools::shift(component.shift_vector, component.triangulation);
  component.offset = offset;
  // Merge the component with the current triangulation
  dealii::distributed::Triangulation<dim> tmp_triangulation(
      _triangulation->get_communicator());
  dealii::GridGenerator::merge_triangulations(
      *_triangulation, component.triangulation, tmp_triangulation);
  _triangulation->clear();
  _triangulation->copy_triangulation(tmp_triangulation);
  // collector_c only needs to be shifted vertically the first time
  if (component.shift_vector[dim - 1] != 0.)
    component.shift_vector[dim - 1] = 0.;
}
}

template <int dim>
unsigned int Geometry<dim>::compute_cell_weight(
    typename dealii::Triangulation<dim, dim>::cell_iterator const &cell)
{
  // Cells in the anode of the cathode have to deal with two physics instead of
  // only one in the collectors and the separator. Each cell starts with a
  // default weight of 1000. This function returns the extra weight on some of
  // the cells.
  unsigned int weight = 0;
  dealii::types::material_id material = cell->material_id();
  std::vector<dealii::types::material_id> &anode = (*_materials)["anode"];
  std::vector<dealii::types::material_id> &cathode = (*_materials)["cathode"];
  if ((std::find(anode.begin(), anode.end(), material) != anode.end()) ||
      (std::find(cathode.begin(), cathode.end(), material) != cathode.end()))
    weight = 1000;

  return weight;
}

template <int dim>
void Geometry<dim>::fill_materials_map(
    std::shared_ptr<boost::property_tree::ptree const> database)
{
  BOOST_ASSERT_MSG(database != nullptr, "The database does not exist.");
  _materials = std::make_shared<std::unordered_map<
      std::string, std::vector<dealii::types::material_id>>>();
  int const n_materials = database->get<int>("materials");
  for (int m = 0; m < n_materials; ++m)
  {
    boost::property_tree::ptree material_database =
        database->get_child("material_" + std::to_string(m));
    std::vector<dealii::types::material_id> const material_ids =
        to_vector<dealii::types::material_id>(
            material_database.get<std::string>("material_id"));
    std::string const material_name =
        material_database.get<std::string>("name");
    _materials->emplace(material_name, material_ids);
  }
}

template <int dim>
void Geometry<dim>::fill_materials_map()
{
  _materials = std::make_shared<std::unordered_map<
      std::string, std::vector<dealii::types::material_id>>>();
  (*_materials)["anode"] = std::vector<dealii::types::material_id>(1, 0);
  (*_materials)["separator"] = std::vector<dealii::types::material_id>(1, 1);
  (*_materials)["cathode"] = std::vector<dealii::types::material_id>(1, 2);
  (*_materials)["collector_anode"] =
      std::vector<dealii::types::material_id>(1, 3);
  (*_materials)["collector_cathode"] =
      std::vector<dealii::types::material_id>(1, 4);
  std::vector<dealii::types::material_id> coll_material_ids(2);
  coll_material_ids[0] = 3;
  coll_material_ids[1] = 4;
  (*_materials)["collector"] = coll_material_ids;
}

template <int dim>
void Geometry<dim>::convert_geometry_database(
    std::shared_ptr<boost::property_tree::ptree> database)
{
  double const cm_to_m = 0.01;
  double const cm2_to_m2 = 0.0001;
  // TODO for now assume that the two collectors have the same dimensions.
  double collector_thickness =
      database->get<double>("anode_collector_thickness") * cm_to_m;
  double anode_thickness =
      database->get<double>("anode_electrode_thickness") * cm_to_m;
  double separator_thickness =
      database->get<double>("separator_thickness") * cm_to_m;
  double cathode_thickness =
      database->get<double>("cathode_electrode_thickness") * cm_to_m;
  double tab_height = database->get<double>("tab_height") * cm_to_m;
  double geometric_area = database->get<double>("geometric_area") * cm2_to_m2;

  if (collector_thickness !=
      database->get<double>("cathode_collector_thickness") * cm_to_m)
    throw std::runtime_error("Both collectors must have the same thickness.");

  if (dim == 2)
  {
    std::string const area_str = std::to_string(geometric_area);
    std::string collector_dim(std::to_string(collector_thickness) + "," +
                              std::to_string(geometric_area + tab_height));
    std::string anode_dim(std::to_string(anode_thickness) + "," + area_str);
    std::string separator_dim(std::to_string(separator_thickness) + "," +
                              area_str);
    std::string cathode_dim(std::to_string(cathode_thickness) + "," + area_str);
    database->put("collector.dimensions", collector_dim);
    database->put("anode.dimensions", anode_dim);
    database->put("separator.dimensions", separator_dim);
    database->put("cathode.dimensions", cathode_dim);
  }
  else
  {
    std::string const area_str = std::to_string(geometric_area);
    std::string collector_dim(std::to_string(collector_thickness) + ",1," +
                              std::to_string(geometric_area + tab_height));
    std::string anode_dim(std::to_string(anode_thickness) + "," + area_str);
    std::string separator_dim(std::to_string(separator_thickness) + ",1," +
                              area_str);
    std::string cathode_dim(std::to_string(cathode_thickness) + ",1," +
                            area_str);
    database->put("collector.dimensions", collector_dim);
    database->put("anode.dimensions", anode_dim);
    database->put("separator.dimensions", separator_dim);
    database->put("cathode.dimensions", cathode_dim);
  }
}

template <int dim>
Geometry<dim>::Geometry(std::shared_ptr<boost::property_tree::ptree> database,
                        boost::mpi::communicator mpi_communicator)
    : _communicator(mpi_communicator),
      _anode_boundary_id(type::invalid_boundary_id),
      _cathode_boundary_id(type::invalid_boundary_id), _triangulation(nullptr),
      _materials(nullptr)
{
  _triangulation = std::make_shared<dealii::distributed::Triangulation<dim>>(
      mpi_communicator);
  std::string mesh_type = database->get<std::string>("type");
  if (mesh_type.compare("file") == 0)
  {
    std::string mesh_file = database->get<std::string>("mesh_file");
    dealii::GridIn<dim> mesh_reader;
    mesh_reader.attach_triangulation(*_triangulation);
    std::fstream fin;
    fin.open(mesh_file.c_str(), std::fstream::in);
    std::string const file_extension =
        mesh_file.substr(mesh_file.find_last_of(".") + 1);
    fill_materials_map(database);
    boost::property_tree::ptree boundary_database =
        database->get_child("boundary_values");
    _anode_boundary_id =
        boundary_database.get<dealii::types::boundary_id>("anode_boundary_id");
    _cathode_boundary_id = boundary_database.get<dealii::types::boundary_id>(
        "cathode_boundary_id");
    if (file_extension.compare("ucd") == 0)
    {
      mesh_reader.read_ucd(fin);
    }
    else if (file_extension.compare("inp") == 0)
    {
      mesh_reader.read_abaqus(fin);
    }
    else
    {
      throw std::runtime_error("Bad mesh file extension ." + file_extension +
                               " in mesh file " + mesh_file);
    }
    fin.close();
  }
  else
  {
    fill_materials_map();
    convert_geometry_database(database);

    // If the mesh type is supercapacitor, we provide a default mesh
    if (mesh_type.compare("supercapacitor") == 0)
    {
      if (dim == 2)
      {
        std::string collector_div("1,6");
        std::string anode_div("10,5");
        std::string separator_div("5,5");
        std::string cathode_div("10,5");
        database->put("collector.divisions", collector_div);
        database->put("anode.divisions", anode_div);
        database->put("separator.divisions", separator_div);
        database->put("cathode.divisions", cathode_div);
      }
      else
      {
        std::string collector_div("3,3,3");
        std::string anode_div("5,5,2");
        std::string separator_div("4,4,2");
        std::string cathode_div("5,5,2");
        database->put("collector.divisions", collector_div);
        database->put("anode.divisions", anode_div);
        database->put("separator.divisions", separator_div);
        database->put("cathode.divisions", cathode_div);
      }

      database->put("n_repetitions", 0);
      database->put("n_refinements", 1);
    }
    mesh_generator(*database);
  }

  // We need to do load balancing because cells in the collectors and the
  // separator don't have both physics.
  _triangulation->signals.cell_weight.connect(std::bind(
      &Geometry<dim>::compute_cell_weight, this, std::placeholders::_1));
  _triangulation->repartition();
}

template <int dim>
Geometry<dim>::Geometry(
    std::shared_ptr<boost::property_tree::ptree const> database,
    std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation)
    : _communicator(boost::mpi::communicator(triangulation->get_communicator(),
                                             boost::mpi::comm_duplicate)),
      _anode_boundary_id(type::invalid_boundary_id),
      _cathode_boundary_id(type::invalid_boundary_id),
      _triangulation(triangulation), _materials(nullptr)
{
  fill_materials_map(database);
}

template <int dim>
void Geometry<dim>::set_boundary_ids(double const collector_top,
                                     double const collector_bottom)
{
  // TODO the code below can be cleaned up when using the next version of
  // deal.II
  //(current is 8.4). Only the for loops will be left without the if. Instead we
  // can use dealii::IteratorFilters::AtBoundary().
  typedef dealii::FilteredIterator<
      typename dealii::Triangulation<dim>::active_cell_iterator> FI;

  // Set the anode boundary id
  bool const locally_owned = true;
  std::set<dealii::types::material_id> collector_anode(
      (*_materials)["collector_anode"].begin(),
      (*_materials)["collector_anode"].end());
  FI anode_cell(dealii::IteratorFilters::MaterialIdEqualTo(collector_anode,
                                                           locally_owned)),
      anode_end_cell(dealii::IteratorFilters::MaterialIdEqualTo(collector_anode,
                                                                locally_owned),
                     _triangulation->end());
  anode_cell.set_to_next_positive(_triangulation->begin_active());
  double const eps = 1e-6;
  unsigned int boundary_id_set = 0;
  for (; anode_cell < anode_end_cell; ++anode_cell)
    if (anode_cell->at_boundary())
      for (unsigned int i = 0; i < dealii::GeometryInfo<dim>::faces_per_cell;
           ++i)
        // Check that the face is on the top of the collector
        if (std::abs(anode_cell->face(i)->center()[dim - 1] - collector_top) <
            eps * anode_cell->measure())
        {
          anode_cell->face(i)->set_boundary_id(_anode_boundary_id);
          boundary_id_set = 1;
        }
  boundary_id_set = dealii::Utilities::MPI::max(boundary_id_set, _communicator);
  BOOST_ASSERT_MSG(boundary_id_set == 1, "Anode boundary id no set.");

  // Set the cathode boundary id
  std::set<dealii::types::material_id> collector_cathode(
      (*_materials)["collector_cathode"].begin(),
      (*_materials)["collector_cathode"].end());
  FI cathode_cell(dealii::IteratorFilters::MaterialIdEqualTo(collector_cathode,
                                                             locally_owned)),
      cathode_end_cell(dealii::IteratorFilters::MaterialIdEqualTo(
                           collector_cathode, locally_owned),
                       _triangulation->end());
  cathode_cell.set_to_next_positive(_triangulation->begin_active());
  boundary_id_set = 0;
  for (; cathode_cell < cathode_end_cell; ++cathode_cell)
    if (cathode_cell->at_boundary())
      for (unsigned int i = 0; i < dealii::GeometryInfo<dim>::faces_per_cell;
           ++i)
        // Check that the face is on the bottom of the collector
        if (std::abs(cathode_cell->face(i)->center()[dim - 1] -
                     collector_bottom) < eps * cathode_cell->measure())
        {
          cathode_cell->face(i)->set_boundary_id(_cathode_boundary_id);
          boundary_id_set = 1;
        }
  boundary_id_set = dealii::Utilities::MPI::max(boundary_id_set, _communicator);
  BOOST_ASSERT_MSG(boundary_id_set == 1, "Cathode boundary id no set.");
}

template <int dim>
void Geometry<dim>::mesh_generator(boost::property_tree::ptree const &database)
{
  // Read the data needed for the collectors
  internal::Component<dim> collector_a(_communicator);
  boost::property_tree::ptree collector_database =
      database.get_child("collector");
  internal::read_component_database(collector_database, collector_a);
  internal::Component<dim> collector_c(collector_a.box_dimensions,
                                       collector_a.repetitions, _communicator);

  // Read the data needed fot the anode
  internal::Component<dim> anode(_communicator);
  boost::property_tree::ptree anode_database = database.get_child("anode");
  internal::read_component_database(anode_database, anode);

  // Read the data needed fot the cathode
  internal::Component<dim> cathode(_communicator);
  boost::property_tree::ptree cathode_database = database.get_child("cathode");
  internal::read_component_database(cathode_database, cathode);

  // Read the data needed fot the separator
  internal::Component<dim> separator(_communicator);
  boost::property_tree::ptree separator_database =
      database.get_child("separator");
  internal::read_component_database(separator_database, separator);

  // For now, we assume that the user does not create hanging nodes with the
  // repetitions
  // Create the triangulation for the anode.
  dealii::GridGenerator::subdivided_hyper_rectangle(
      anode.triangulation, anode.repetitions, anode.box_dimensions[0],
      anode.box_dimensions[1]);
  for (auto cell : anode.triangulation.cell_iterators())
    cell->set_material_id((*_materials)["anode"][0]);
  // Create the triangulation for the cathode.
  dealii::GridGenerator::subdivided_hyper_rectangle(
      cathode.triangulation, cathode.repetitions, cathode.box_dimensions[0],
      cathode.box_dimensions[1]);
  for (auto cell : cathode.triangulation.cell_iterators())
    cell->set_material_id((*_materials)["cathode"][0]);
  // Create the triangulation for the seperator.
  dealii::GridGenerator::subdivided_hyper_rectangle(
      separator.triangulation, separator.repetitions,
      separator.box_dimensions[0], separator.box_dimensions[1]);
  for (auto cell : separator.triangulation.cell_iterators())
    cell->set_material_id((*_materials)["separator"][0]);

  // Create the triangulation for first collector.
  double const anode_dim = anode.box_dimensions[1][dim - 1];
  double const collector_dim = collector_a.box_dimensions[1][dim - 1];
  double const delta_collector =
      collector_dim / collector_a.repetitions[dim - 1];
  dealii::GridGenerator::subdivided_hyper_rectangle(
      collector_a.triangulation, collector_a.repetitions,
      collector_a.box_dimensions[0], collector_a.box_dimensions[1]);
  for (auto cell : collector_a.triangulation.cell_iterators())
    cell->set_material_id((*_materials)["collector_anode"][0]);
  double const scale_factor_a = anode_dim / (collector_dim - delta_collector);
  std::function<dealii::Point<dim>(dealii::Point<dim> const &)> transform_a =
      std::bind(&internal::transform_coll_a<dim>, std::placeholders::_1,
                scale_factor_a, collector_a.box_dimensions[1][dim - 1]);
  dealii::GridTools::transform(transform_a, collector_a.triangulation);

  // Create the triangulation for the second collector. For now, we assume that
  // collector_a and collector_c have the same mesh.
  dealii::GridGenerator::subdivided_hyper_rectangle(
      collector_c.triangulation, collector_c.repetitions,
      collector_c.box_dimensions[0], collector_c.box_dimensions[1]);
  for (auto cell : collector_c.triangulation.cell_iterators())
    cell->set_material_id((*_materials)["collector_cathode"][0]);
  double const scale_factor_c = anode_dim / (collector_dim - delta_collector);
  std::function<dealii::Point<dim>(dealii::Point<dim> const &)> transform_c =
      std::bind(&internal::transform_coll_c<dim>, std::placeholders::_1,
                scale_factor_c, collector_c.box_dimensions[1][dim - 1],
                collector_dim - anode_dim - scale_factor_c * delta_collector);
  dealii::GridTools::transform(transform_c, collector_c.triangulation);
  collector_c.shift_vector[dim - 1] = -(collector_dim - anode_dim);

  std::array<internal::Component<dim> *, 8> components = {
      &anode,   &separator, &cathode, &collector_c,
      &cathode, &separator, &anode,   &collector_a};
  _triangulation->clear();
  _triangulation->copy_triangulation(collector_a.triangulation);
  double offset = collector_a.box_dimensions[1][0];
  unsigned int pos = 0;
  unsigned int const n_repetitions =
      database.get<unsigned int>("n_repetitions", 1);
  for (unsigned int i = 0; i <= n_repetitions; ++i)
  {
    for (unsigned int j = 0; j < 4; ++j)
    {
      internal::merge_components(*components[pos], offset, _triangulation);
      offset += components[pos]->box_dimensions[1][0];
      ++pos;
      pos = pos % 8;
    }
  }

  // Apply global refinement
  unsigned int const n_refinements =
      database.get<unsigned int>("n_refinements", 0);
  _triangulation->refine_global(n_refinements);

  // Apply boundary conditions. This needs to be done after the merging because
  // the merging loses the boundary id
  _anode_boundary_id = 1;
  _cathode_boundary_id = 2;
  set_boundary_ids(collector_a.box_dimensions[1][dim - 1],
                   -(collector_dim - anode_dim));
}

} // end namespace cap
