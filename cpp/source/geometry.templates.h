#include <cap/geometry.h>
#include <cap/utils.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/geometry_info.h>
#include <fstream>
#include <tuple>

namespace cap
{

template <int dim>
std::tuple<std::pair<dealii::Point<dim>, dealii::Point<dim>>,
           std::pair<dealii::Point<dim>, dealii::Point<dim>>,
           std::pair<dealii::Point<dim>, dealii::Point<dim>>,
           std::pair<dealii::Point<dim>, dealii::Point<dim>>,
           std::pair<dealii::Point<dim>, dealii::Point<dim>>,
           std::pair<dealii::Point<dim>, dealii::Point<dim>>,
           std::pair<dealii::Point<dim>, dealii::Point<dim>>>
foo(std::shared_ptr<boost::property_tree::ptree const> const &database);

template <>
std::tuple<std::pair<dealii::Point<2>, dealii::Point<2>>,
           std::pair<dealii::Point<2>, dealii::Point<2>>,
           std::pair<dealii::Point<2>, dealii::Point<2>>,
           std::pair<dealii::Point<2>, dealii::Point<2>>,
           std::pair<dealii::Point<2>, dealii::Point<2>>,
           std::pair<dealii::Point<2>, dealii::Point<2>>,
           std::pair<dealii::Point<2>, dealii::Point<2>>>
foo<2>(std::shared_ptr<boost::property_tree::ptree const> const &database)
{
  // clang-format off
  double const anode_collector_width   = database->get<double>("anode_collector_width");
  double const anode_electrode_width   = database->get<double>("anode_electrode_width");
  double const separator_width         = database->get<double>("separator_width");
  double const cathode_electrode_width = database->get<double>("cathode_electrode_width");
  double const cathode_collector_width = database->get<double>("cathode_collector_width");
  double const sandwich_height         = database->get<double>("sandwich_height");
  double const tab_height              = database->get<double>("tab_height");
  // clang-format on

  std::pair<dealii::Point<2>, dealii::Point<2>> anode_tab_bbox;
  std::pair<dealii::Point<2>, dealii::Point<2>> anode_collector_bbox;
  std::pair<dealii::Point<2>, dealii::Point<2>> anode_electrode_bbox;
  std::pair<dealii::Point<2>, dealii::Point<2>> separator_bbox;
  std::pair<dealii::Point<2>, dealii::Point<2>> cathode_electrode_bbox;
  std::pair<dealii::Point<2>, dealii::Point<2>> cathode_collector_bbox;
  std::pair<dealii::Point<2>, dealii::Point<2>> cathode_tab_bbox;

  anode_tab_bbox =
      std::make_pair(dealii::Point<2>(-anode_collector_width, sandwich_height),
                     dealii::Point<2>(0.0, sandwich_height + tab_height));
  anode_collector_bbox =
      std::make_pair(dealii::Point<2>(-anode_collector_width, 0.0),
                     dealii::Point<2>(0.0, sandwich_height));
  anode_electrode_bbox =
      std::make_pair(dealii::Point<2>(0.0, 0.0),
                     dealii::Point<2>(anode_electrode_width, sandwich_height));
  separator_bbox =
      std::make_pair(dealii::Point<2>(anode_electrode_width, 0.0),
                     dealii::Point<2>(anode_electrode_width + separator_width,
                                      sandwich_height));
  cathode_electrode_bbox = std::make_pair(
      dealii::Point<2>(anode_electrode_width + separator_width, 0.0),
      dealii::Point<2>(anode_electrode_width + separator_width +
                           cathode_electrode_width,
                       sandwich_height));
  cathode_collector_bbox = std::make_pair(
      dealii::Point<2>(anode_electrode_width + separator_width +
                           cathode_electrode_width,
                       0.0),
      dealii::Point<2>(anode_electrode_width + separator_width +
                           cathode_electrode_width + cathode_collector_width,
                       sandwich_height));
  cathode_tab_bbox = std::make_pair(
      dealii::Point<2>(anode_electrode_width + separator_width +
                           cathode_electrode_width,
                       -tab_height),
      dealii::Point<2>(anode_electrode_width + separator_width +
                           cathode_electrode_width + cathode_collector_width,
                       0.0));
  return std::make_tuple(anode_tab_bbox, anode_collector_bbox,
                         anode_electrode_bbox, separator_bbox,
                         cathode_electrode_bbox, cathode_collector_bbox,
                         cathode_tab_bbox);
}

template <>
std::tuple<std::pair<dealii::Point<3>, dealii::Point<3>>,
           std::pair<dealii::Point<3>, dealii::Point<3>>,
           std::pair<dealii::Point<3>, dealii::Point<3>>,
           std::pair<dealii::Point<3>, dealii::Point<3>>,
           std::pair<dealii::Point<3>, dealii::Point<3>>,
           std::pair<dealii::Point<3>, dealii::Point<3>>,
           std::pair<dealii::Point<3>, dealii::Point<3>>>
foo<3>(std::shared_ptr<boost::property_tree::ptree const> const &database)
{
  // clang-format off
  double const anode_collector_width   = database->get<double>("anode_collector_width");
  double const anode_electrode_width   = database->get<double>("anode_electrode_width");
  double const separator_width         = database->get<double>("separator_width");
  double const cathode_electrode_width = database->get<double>("cathode_electrode_width");
  double const cathode_collector_width = database->get<double>("cathode_collector_width");
  double const sandwich_height         = database->get<double>("sandwich_height");
  double const sandwich_depth          = database->get<double>("sandwich_depth");
  double const tab_height              = database->get<double>("tab_height");
  // clang-format on

  std::pair<dealii::Point<3>, dealii::Point<3>> anode_tab_bbox;
  std::pair<dealii::Point<3>, dealii::Point<3>> anode_collector_bbox;
  std::pair<dealii::Point<3>, dealii::Point<3>> anode_electrode_bbox;
  std::pair<dealii::Point<3>, dealii::Point<3>> separator_bbox;
  std::pair<dealii::Point<3>, dealii::Point<3>> cathode_electrode_bbox;
  std::pair<dealii::Point<3>, dealii::Point<3>> cathode_collector_bbox;
  std::pair<dealii::Point<3>, dealii::Point<3>> cathode_tab_bbox;

  anode_tab_bbox = std::make_pair(
      dealii::Point<3>(-anode_collector_width, sandwich_height, 0.0),
      dealii::Point<3>(0.0, sandwich_height + tab_height, sandwich_depth));
  anode_collector_bbox =
      std::make_pair(dealii::Point<3>(-anode_collector_width, 0.0, 0.0),
                     dealii::Point<3>(0.0, sandwich_height, sandwich_depth));
  anode_electrode_bbox = std::make_pair(
      dealii::Point<3>(0.0, 0.0, 0.0),
      dealii::Point<3>(anode_electrode_width, sandwich_height, sandwich_depth));
  separator_bbox =
      std::make_pair(dealii::Point<3>(anode_electrode_width, 0.0, 0.0),
                     dealii::Point<3>(anode_electrode_width + separator_width,
                                      sandwich_height, sandwich_depth));
  cathode_electrode_bbox = std::make_pair(
      dealii::Point<3>(anode_electrode_width + separator_width, 0.0, 0.0),
      dealii::Point<3>(anode_electrode_width + separator_width +
                           cathode_electrode_width,
                       sandwich_height, sandwich_depth));
  cathode_collector_bbox = std::make_pair(
      dealii::Point<3>(anode_electrode_width + separator_width +
                           cathode_electrode_width,
                       0.0, 0.0),
      dealii::Point<3>(anode_electrode_width + separator_width +
                           cathode_electrode_width + cathode_collector_width,
                       sandwich_height, sandwich_depth));
  cathode_tab_bbox = std::make_pair(
      dealii::Point<3>(anode_electrode_width + separator_width +
                           cathode_electrode_width,
                       -tab_height, 0.0),
      dealii::Point<3>(anode_electrode_width + separator_width +
                           cathode_electrode_width + cathode_collector_width,
                       0.0, sandwich_depth));
  return std::make_tuple(anode_tab_bbox, anode_collector_bbox,
                         anode_electrode_bbox, separator_bbox,
                         cathode_electrode_bbox, cathode_collector_bbox,
                         cathode_tab_bbox);
}

template <int dim>
Geometry<dim>::Geometry(
    std::shared_ptr<boost::property_tree::ptree const> const &database)
{
  this->triangulation         = std::make_shared<dealii::Triangulation<dim>>();
  std::string const mesh_file = database->get<std::string>("mesh_file");
  dealii::GridIn<dim> mesh_reader;
  mesh_reader.attach_triangulation(*(this->triangulation));
  std::fstream fin;
  fin.open(mesh_file.c_str(), std::fstream::in);
  std::string const file_extension =
      mesh_file.substr(mesh_file.find_last_of(".") + 1);
  if (file_extension.compare("ucd") == 0)
  {
    mesh_reader.read_ucd(fin);
  }
  else
  {
    throw std::runtime_error("Bad mesh file extension ." + file_extension +
                             " in mesh file " + mesh_file);
  }
  fin.close();

  this->materials = std::make_shared<std::unordered_map<
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
    (*this->materials).emplace(material_name, material_ids);
  }
}

template <int dim>
SuperCapacitorGeometry<dim>::SuperCapacitorGeometry(
    std::shared_ptr<boost::property_tree::ptree const> const &database)
    : Geometry<dim>(database)
{
  // clang-format off
  this->separator_material_id         = database->get<dealii::types::material_id>("separator_material_id");
  this->anode_electrode_material_id   = database->get<dealii::types::material_id>("anode_electrode_material_id");
  this->anode_collector_material_id   = database->get<dealii::types::material_id>("anode_collector_material_id");
  this->cathode_electrode_material_id =database->get<dealii::types::material_id>("cathode_electrode_material_id");
  this->cathode_collector_material_id = database->get<dealii::types::material_id>("cathode_collector_material_id");
  // clang-format on

  std::unordered_map<dealii::types::material_id,
                     std::pair<dealii::Point<dim>, dealii::Point<dim>>> bboxes;
  auto initialize_bbox = [](void)
  {
    double const max = std::numeric_limits<double>::max();
    double const min = std::numeric_limits<double>::lowest();
    if (dim > 2)
      return std::make_pair(dealii::Point<dim>(max, max, max),
                            dealii::Point<dim>(min, min, min));
    else
      return std::make_pair(dealii::Point<dim>(max, max),
                            dealii::Point<dim>(min, min));
  };
  for (dealii::types::material_id material_id : {
           this->separator_material_id, this->anode_electrode_material_id,
           this->anode_collector_material_id,
           this->cathode_electrode_material_id,
           this->cathode_collector_material_id,
       })
    bboxes.emplace(material_id, initialize_bbox());
  auto add_vertex_to_bbox =
      [](dealii::Point<dim> const &p,
         std::pair<dealii::Point<dim>, dealii::Point<dim>> &bbox)
  {
    for (int d = 0; d < dim; ++d)
    {
      if (p[d] < bbox.first[d])
        bbox.first[d] = p[d];
      if (p[d] > bbox.second[d])
        bbox.second[d] = p[d];
    }
  };
  int const vertices_per_cell = dealii::GeometryInfo<dim>::vertices_per_cell;
  auto cell                   = (*this->triangulation).begin_active();
  auto end_cell = (*this->triangulation).end();
  for (; cell != end_cell; ++cell)
  {
    auto material_id = cell->material_id();
    for (int vertex = 0; vertex < vertices_per_cell; ++vertex)
      add_vertex_to_bbox(cell->vertex(vertex), bboxes[material_id]);
  }

  this->anode_tab_bbox         = bboxes[this->anode_collector_material_id];
  this->anode_collector_bbox   = bboxes[this->anode_collector_material_id];
  this->anode_electrode_bbox   = bboxes[this->anode_electrode_material_id];
  this->separator_bbox         = bboxes[this->separator_material_id];
  this->cathode_electrode_bbox = bboxes[this->cathode_electrode_material_id];
  this->cathode_collector_bbox = bboxes[this->cathode_collector_material_id];
  this->cathode_tab_bbox = bboxes[this->cathode_collector_material_id];
  AssertThrow((this->anode_electrode_bbox).second[0] ==
                  (this->separator_bbox).first[0],
              dealii::StandardExceptions::ExcMessage("duh"));
  AssertThrow((this->anode_electrode_bbox).first[1] ==
                  (this->separator_bbox).first[1],
              dealii::StandardExceptions::ExcMessage("duh"));
  AssertThrow((this->separator_bbox).second[0] ==
                  (this->cathode_electrode_bbox).first[0],
              dealii::StandardExceptions::ExcMessage("duh"));
  AssertThrow((this->separator_bbox).first[1] ==
                  (this->cathode_electrode_bbox).first[1],
              dealii::StandardExceptions::ExcMessage("duh"));
  (this->anode_tab_bbox).first[1] = (this->anode_electrode_bbox).second[1];
  (this->anode_collector_bbox).second[1] =
      (this->anode_electrode_bbox).second[1];
  AssertThrow((this->anode_collector_bbox).second[0] ==
                  (this->anode_electrode_bbox).first[0],
              dealii::StandardExceptions::ExcMessage("duh"));
  AssertThrow((this->anode_collector_bbox).first[1] ==
                  (this->anode_electrode_bbox).first[1],
              dealii::StandardExceptions::ExcMessage("duh"));
  (this->cathode_collector_bbox).first[1] =
      (this->cathode_electrode_bbox).first[1];
  (this->cathode_tab_bbox).second[1] = (this->cathode_electrode_bbox).first[1];
  AssertThrow((this->cathode_electrode_bbox).second[0] ==
                  (this->cathode_collector_bbox).first[0],
              dealii::StandardExceptions::ExcMessage("duh"));
  AssertThrow((this->cathode_electrode_bbox).first[1] ==
                  (this->cathode_collector_bbox).first[1],
              dealii::StandardExceptions::ExcMessage("duh"));

  this->reset(database);
}

template <int dim>
void SuperCapacitorGeometry<dim>::reset(
    std::shared_ptr<boost::property_tree::ptree const> const &database)
{
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_anode_tab_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_anode_collector_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_anode_electrode_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_separator_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_cathode_electrode_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_cathode_collector_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> new_cathode_tab_bbox;
  std::tie(new_anode_tab_bbox, new_anode_collector_bbox,
           new_anode_electrode_bbox, new_separator_bbox,
           new_cathode_electrode_bbox, new_cathode_collector_bbox,
           new_cathode_tab_bbox) = foo<dim>(database);

  auto transform_point =
      [](dealii::Point<dim> const &p_0, dealii::Point<dim> const &p_1,
         dealii::Point<dim> const &p, dealii::Point<dim> const &q_0,
         dealii::Point<dim> const &q_1, dealii::Point<dim> &q)
  {
    for (int d = 0; d < dim; ++d)
    {
      q[d] = (q_1[d] - q_0[d]) / (p_1[d] - p_0[d]) * (p[d] - p_0[d]) + q_0[d];
    } // end for d
  };
  std::vector<bool> vertex_visited((*this->triangulation).n_vertices(), false);
  typedef typename dealii::Triangulation<dim, spacedim>::active_cell_iterator
      active_cell_iterator;
  auto move_vertices_in_cell = [&transform_point, &vertex_visited](
      active_cell_iterator &cell,
      std::pair<dealii::Point<dim>, dealii::Point<dim>> const &old_bbox,
      std::pair<dealii::Point<dim>, dealii::Point<dim>> const &new_bbox)
  {
    for (unsigned int vertex = 0;
         vertex < dealii::GeometryInfo<dim>::vertices_per_cell; ++vertex)
    {
      if (!vertex_visited[cell->vertex_index(vertex)])
      {
        transform_point(old_bbox.first, old_bbox.second, cell->vertex(vertex),
                        new_bbox.first, new_bbox.second, cell->vertex(vertex));
        vertex_visited[cell->vertex_index(vertex)] = true;
      } // end if vertex
    }   // end for vertex
  };
  auto point_in_bbox =
      [](dealii::Point<dim> const &p,
         std::pair<dealii::Point<dim>, dealii::Point<dim>> const &bbox)
  {
    double const tol = 1.0e-10 * (bbox.second - bbox.first).norm();
    for (int d = 0; d < dim; ++d)
      if ((bbox.first[d] - tol > p[d]) || (p[d] > bbox.second[d] + tol))
        return false;
    return true;
  };

  for (auto cell = (*this->triangulation).begin_active();
       cell != (*this->triangulation).end(); ++cell)
    if (point_in_bbox(cell->center(), this->anode_tab_bbox) ||
        point_in_bbox(cell->center(), this->cathode_tab_bbox))
      cell->set_user_flag();

  // clang-format off
  std::pair<dealii::Point<dim>, dealii::Point<dim>> const &old_anode_collector_bbox   = this->anode_collector_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> const &old_anode_tab_bbox         = this->anode_tab_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> const &old_cathode_collector_bbox = this->cathode_collector_bbox;
  std::pair<dealii::Point<dim>, dealii::Point<dim>> const &old_cathode_tab_bbox       = this->cathode_tab_bbox;
  // clang-format on

  std::map<dealii::types::material_id,
           std::function<void(active_cell_iterator &)>> cell_transform;
  cell_transform[this->anode_collector_material_id] =
      [&move_vertices_in_cell, &old_anode_tab_bbox, &new_anode_tab_bbox,
       &old_anode_collector_bbox,
       &new_anode_collector_bbox](active_cell_iterator &cell)
  {
    if (cell->user_flag_set())
      move_vertices_in_cell(cell, old_anode_tab_bbox, new_anode_tab_bbox);
    else
      move_vertices_in_cell(cell, old_anode_collector_bbox,
                            new_anode_collector_bbox);
  };
  cell_transform[this->anode_electrode_material_id] =
      std::bind(move_vertices_in_cell, std::placeholders::_1,
                this->anode_electrode_bbox, new_anode_electrode_bbox);
  cell_transform[this->separator_material_id] =
      std::bind(move_vertices_in_cell, std::placeholders::_1,
                this->separator_bbox, new_separator_bbox);
  cell_transform[this->cathode_electrode_material_id] =
      std::bind(move_vertices_in_cell, std::placeholders::_1,
                this->cathode_electrode_bbox, new_cathode_electrode_bbox);
  cell_transform[this->cathode_collector_material_id] =
      [&move_vertices_in_cell, &old_cathode_tab_bbox, &new_cathode_tab_bbox,
       &old_cathode_collector_bbox,
       &new_cathode_collector_bbox](active_cell_iterator &cell)
  {
    if (cell->user_flag_set())
      move_vertices_in_cell(cell, old_cathode_tab_bbox, new_cathode_tab_bbox);
    else
      move_vertices_in_cell(cell, old_cathode_collector_bbox,
                            new_cathode_collector_bbox);
  };

  typename dealii::Triangulation<dim>::active_cell_iterator cell =
      (*this->triangulation).begin_active();
  typename dealii::Triangulation<dim>::active_cell_iterator end_cell =
      (*this->triangulation).end();
  for (; cell != end_cell; ++cell)
  {
    bool found_it = false;
    for (typename std::map<
             dealii::types::material_id,
             std::function<void(active_cell_iterator & cell)>>::iterator it =
             cell_transform.begin();
         it != cell_transform.end(); ++it)
    {
      if (it->first == cell->material_id())
      {
        it->second(cell);
        found_it = true;
      }
    }
    if (!found_it)
      throw std::runtime_error("Error while moving the vertices");
  } // end for cell
  this->anode_tab_bbox         = new_anode_tab_bbox;
  this->anode_collector_bbox   = new_anode_collector_bbox;
  this->anode_electrode_bbox   = new_anode_electrode_bbox;
  this->separator_bbox         = new_separator_bbox;
  this->cathode_collector_bbox = new_cathode_collector_bbox;
  this->cathode_electrode_bbox = new_cathode_electrode_bbox;
  this->cathode_tab_bbox       = new_cathode_tab_bbox;

  (*this->triangulation).clear_user_flags();
}

} // end namespace cap
