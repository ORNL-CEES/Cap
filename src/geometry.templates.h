#include <cap/geometry.h>
#include <deal.II/grid/grid_in.h>
#include <fstream>

namespace cap {

template <int dim>
Geometry<dim>::
Geometry(std::shared_ptr<boost::property_tree::ptree const> const & database)
{
    this->triangulation = std::make_shared<dealii::Triangulation<dim> >();
    std::string const mesh_file = database->get<std::string>("mesh_file");
    dealii::GridIn<dim> mesh_reader;
    mesh_reader.attach_triangulation(*(this->triangulation));
    std::fstream fin;
    fin.open(mesh_file.c_str(), std::fstream::in);
    std::string const file_extension = mesh_file.substr(mesh_file.find_last_of(".")+1);
    if (file_extension.compare("ucd") == 0) {
        mesh_reader.read_ucd(fin);
    } else {
        throw std::runtime_error("Bad mesh file extension ."+file_extension+" in mesh file "+mesh_file);
    }
    fin.close();
}



template <int dim>
bool point_in_bbox(dealii::Point<dim> const & p, 
    std::pair<dealii::Point<dim>,dealii::Point<dim> > const & bbox) 
{
    double const tol = 1.0e-10 * (bbox.second - bbox.first).norm();
    for (int d = 0; d < dim; ++d)
        if ((bbox.first[d] - tol > p[d]) || (p[d] > bbox.second[d] + tol))
            return false;
    return true;
}



template <int dim>
SuperCapacitorGeometry<dim>::
SuperCapacitorGeometry(std::shared_ptr<boost::property_tree::ptree const> const & database)
: Geometry<dim>(database)
{
    this->separator_material_id          = database->get<dealii::types::material_id>("separator_material_id"        );
    this->anode_electrode_material_id    = database->get<dealii::types::material_id>("anode_electrode_material_id"  );
    this->anode_collector_material_id    = database->get<dealii::types::material_id>("anode_collector_material_id"  );
    this->cathode_electrode_material_id  = database->get<dealii::types::material_id>("cathode_electrode_material_id");
    this->cathode_collector_material_id  = database->get<dealii::types::material_id>("cathode_collector_material_id");
    double const electrode_width         = database->get<double                    >("electrode_width"              );
    double const separator_width         = database->get<double                    >("separator_width"              );
    double const collector_width         = database->get<double                    >("collector_width"              );
    double const sandwich_height         = database->get<double                    >("sandwich_height"              );
    double const tab_height              = database->get<double                    >("tab_height"                   ); 


    this->anode_tab_bbox = 
        std::make_pair(
            dealii::Point<dim>(-collector_width                                               , sandwich_height           ), 
            dealii::Point<dim>(0.0                                                            , sandwich_height+tab_height)
        );
    this->anode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(-collector_width                                               , 0.0                       ), 
            dealii::Point<dim>(0.0                                                            , sandwich_height           )
        );
    this->anode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(0.0                                                            , 0.0                       ), 
            dealii::Point<dim>(electrode_width                                                , sandwich_height           )
        );
    this->separator_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width                                                , 0.0                       ), 
            dealii::Point<dim>(electrode_width+separator_width                                , sandwich_height           )
        );
    this->cathode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width+separator_width                                , 0.0                       ), 
            dealii::Point<dim>(electrode_width+separator_width+electrode_width                , sandwich_height           )
        );
    this->cathode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width+separator_width+electrode_width                , 0.0                       ), 
            dealii::Point<dim>(electrode_width+separator_width+electrode_width+collector_width, sandwich_height           )
        );
    this->cathode_tab_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width+separator_width+electrode_width                , -tab_height               ), 
            dealii::Point<dim>(electrode_width+separator_width+electrode_width+collector_width, 0.0                       )
        );
        
    std::map<dealii::types::material_id, std::function<bool(dealii::Point<dim> const &)> > point_in_material_id;
    point_in_material_id[this->anode_collector_material_id] =
        std::bind(
            std::logical_or<bool>(),
            std::bind(point_in_bbox<dim>, std::placeholders::_1, anode_collector_bbox),
            std::bind(point_in_bbox<dim>, std::placeholders::_1, anode_tab_bbox)
        );
    point_in_material_id[this->anode_electrode_material_id] =
        std::bind(point_in_bbox<dim>, std::placeholders::_1, anode_electrode_bbox);
    point_in_material_id[this->separator_material_id] =
        std::bind(point_in_bbox<dim>, std::placeholders::_1, separator_bbox);
    point_in_material_id[this->cathode_electrode_material_id] =
        std::bind(point_in_bbox<dim>, std::placeholders::_1, cathode_electrode_bbox);
    point_in_material_id[this->cathode_collector_material_id] =
        std::bind(
            std::logical_or<bool>(),
            std::bind(point_in_bbox<dim>, std::placeholders::_1, cathode_collector_bbox),
            std::bind(point_in_bbox<dim>, std::placeholders::_1, cathode_tab_bbox)
        );

    typename dealii::Triangulation<dim>::active_cell_iterator cell     = (this->triangulation)->begin_active();
    typename dealii::Triangulation<dim>::active_cell_iterator end_cell = (this->triangulation)->end();
    for ( ; cell != end_cell; ++cell) {
        bool found_it = false;
        for (typename std::map<dealii::types::material_id, std::function<bool(dealii::Point<dim> const &)> >::const_iterator it = point_in_material_id.begin();
            it != point_in_material_id.end();
            ++it)
        {
            if (it->second(cell->center())) {
                found_it = true;
                cell->set_material_id(it->first);  
            }
        }
        if (!found_it)
            throw std::runtime_error("Error while setting material ids");
    }
}



template <int dim>
void
SuperCapacitorGeometry<dim>::
reset(std::shared_ptr<boost::property_tree::ptree const> const & database)
{
    double const electrode_width = database->get<double>("electrode_width");
    double const separator_width = database->get<double>("separator_width");
    double const collector_width = database->get<double>("collector_width");
    double const sandwich_height = database->get<double>("sandwich_height");
    double const tab_height      = database->get<double>("tab_height"     ); 

    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_anode_tab_bbox = 
        std::make_pair(
            dealii::Point<dim>(-collector_width                                               , sandwich_height           ), 
            dealii::Point<dim>(0.0                                                            , sandwich_height+tab_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_anode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(-collector_width                                               , 0.0                       ), 
            dealii::Point<dim>(0.0                                                            , sandwich_height           )
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_anode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(0.0                                                            , 0.0                       ), 
            dealii::Point<dim>(electrode_width                                                , sandwich_height           )
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_separator_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width                                                , 0.0                       ), 
            dealii::Point<dim>(electrode_width+separator_width                                , sandwich_height           )
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_cathode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width+separator_width                                , 0.0                       ), 
            dealii::Point<dim>(electrode_width+separator_width+electrode_width                , sandwich_height           )
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_cathode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width+separator_width+electrode_width                , 0.0                       ), 
            dealii::Point<dim>(electrode_width+separator_width+electrode_width+collector_width, sandwich_height           )
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_cathode_tab_bbox =
        std::make_pair(
            dealii::Point<dim>(electrode_width+separator_width+electrode_width                , -tab_height               ), 
            dealii::Point<dim>(electrode_width+separator_width+electrode_width+collector_width, 0.0                       )
        );
    
    auto transform_point =
        [](dealii::Point<dim> const & p_0, dealii::Point<dim> const & p_1, dealii::Point<dim> const & p,
           dealii::Point<dim> const & q_0, dealii::Point<dim> const & q_1, dealii::Point<dim> &       q)
        {
            for (int d = 0; d < dim; ++d) {
                q[d] = (q_1[d] - q_0[d]) / (p_1[d] - p_0[d]) * (p[d] - p_0[d]) + q_0[d];
            } // end for d
        };
    std::vector<bool> vertex_visited((this->triangulation)->n_vertices(), false);
    typedef typename dealii::Triangulation<dim, spacedim>::active_cell_iterator active_cell_iterator;
    auto move_vertices_in_cell =
        [& transform_point, & vertex_visited](active_cell_iterator & cell,
           std::pair<dealii::Point<dim>,dealii::Point<dim> > const & old_bbox,
           std::pair<dealii::Point<dim>,dealii::Point<dim> > const & new_bbox)
        {
            for (unsigned int vertex = 0; vertex < dealii::GeometryInfo<dim>::vertices_per_cell; ++vertex) {
                if (!vertex_visited[cell->vertex_index(vertex)]) {
                    transform_point(old_bbox.first, old_bbox.second, cell->vertex(vertex), new_bbox.first, new_bbox.second, cell->vertex(vertex));
                    vertex_visited[cell->vertex_index(vertex)] = true;
                } // end if vertex 
            } // end for vertex
        };
    auto conditional_move =
        [](active_cell_iterator & cell, 
            std::function<bool(active_cell_iterator const &)> const & cond,
            std::function<void(active_cell_iterator &)> const & foo,
            std::function<void(active_cell_iterator &)> const & bar)
        {
            if (cond(cell))
                foo(cell);
            else
                bar(cell);
        };
    auto cell_in_bbox =
        [](active_cell_iterator const & cell,
            std::pair<dealii::Point<dim>,dealii::Point<dim> > bbox)
//            std::pair<dealii::Point<dim>,dealii::Point<dim> > const & bbox)
        {
            return point_in_bbox(cell->center(), bbox);
        };
    auto always_true = [] (active_cell_iterator const &) { return true; };
    auto do_nothing = [] (active_cell_iterator &) { };

    std::map<dealii::types::material_id, std::function<void(active_cell_iterator &)> >  cell_transform;
    cell_transform[this->anode_collector_material_id] =
        std::bind(
            move_vertices_in_cell, std::placeholders::_1,
            this->anode_collector_bbox,
            new_anode_collector_bbox
        );
//        std::bind(
//            conditional_move,
//            std::placeholders::_1,
//            always_true,
////            std::bind(cell_in_bbox, std::placeholders::_1, this->anode_collector_bbox),
//            do_nothing,
//            do_nothing
////            std::bind(
////                move_vertices_in_cell, std::placeholders::_1,
////                this->anode_collector_bbox,
////                new_anode_collector_bbox
////            ),
////            std::bind(
////                move_vertices_in_cell, std::placeholders::_1,
////                this->anode_tab_bbox,
////                new_anode_tab_bbox
////            )
//        );
    cell_transform[this->anode_electrode_material_id] =
        std::bind(
            move_vertices_in_cell, std::placeholders::_1,
            this->anode_electrode_bbox,
            new_anode_electrode_bbox
        );
    cell_transform[this->separator_material_id] =
        std::bind(
            move_vertices_in_cell, std::placeholders::_1,
            this->separator_bbox,
            new_separator_bbox
        );
    cell_transform[this->cathode_electrode_material_id] =
        std::bind(
            move_vertices_in_cell, std::placeholders::_1,
            this->cathode_electrode_bbox,
            new_cathode_electrode_bbox
        );
    cell_transform[this->cathode_collector_material_id] =
        std::bind(
            move_vertices_in_cell, std::placeholders::_1,
            this->cathode_collector_bbox,
            new_cathode_collector_bbox
        );
            
    typename dealii::Triangulation<dim>::active_cell_iterator cell     = (this->triangulation)->begin_active();
    typename dealii::Triangulation<dim>::active_cell_iterator end_cell = (this->triangulation)->end();
    for ( ; cell != end_cell; ++cell) {
        bool found_it = false;
        for (typename std::map<dealii::types::material_id, std::function<void(active_cell_iterator & cell)> >::iterator it = cell_transform.begin();
            it != cell_transform.end();
            ++it)
        {
            if (it->first == cell->material_id()) {
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
}

} // end namespace cap
