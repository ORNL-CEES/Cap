#include <cap/utils.h>
#include <deal.II/base/types.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <fstream>

namespace cap {
template <int dim>
class Geometry {
public:
    Geometry(std::shared_ptr<boost::property_tree::ptree const> const & database);
    inline dealii::Triangulation<dim> const & get_triangulation() const
        { return this->triangulation; }
    void move_mesh(std::shared_ptr<boost::property_tree::ptree const> const & database);

private:
    dealii::types::material_id separator_material_id        ;
    dealii::types::material_id anode_electrode_material_id  ;
    dealii::types::material_id anode_collector_material_id  ;
    dealii::types::material_id cathode_electrode_material_id;
    dealii::types::material_id cathode_collector_material_id;
    double electrode_width;
    double separator_width;
    double collector_width;
    double sandwich_height;
    double tab_height;
    void read_mesh(std::string const & mesh_file);
    void write_mesh(std::string const & mesh_file) const;
    dealii::Triangulation<dim> triangulation;

    std::pair<dealii::Point<dim>,dealii::Point<dim> > anode_tab_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > anode_collector_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > anode_electrode_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > separator_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > cathode_electrode_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > cathode_collector_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > cathode_tab_bbox;
    static int const spacedim = dim;
};

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

//template <int dim>
//bool cell_in_bbox(typename dealii::Triangulation<dim>::active_cell_iterator const & cell,
//    std::pair<dealii::Point<dim>,dealii::Point<dim> > const & bbox)
//{
//    return point_in_bbox<dim>(cell->center(), bbox);
//}

template <int dim>
Geometry<dim>::
Geometry(std::shared_ptr<boost::property_tree::ptree const> const & database)
{
    this->separator_material_id          = database->get<dealii::types::material_id>("separator_material_id"        );
    this->anode_electrode_material_id    = database->get<dealii::types::material_id>("anode_electrode_material_id"  );
    this->anode_collector_material_id    = database->get<dealii::types::material_id>("anode_collector_material_id"  );
    this->cathode_electrode_material_id  = database->get<dealii::types::material_id>("cathode_electrode_material_id");
    this->cathode_collector_material_id  = database->get<dealii::types::material_id>("cathode_collector_material_id");
    this->electrode_width = database->get<double>("electrode_width");
    this->separator_width = database->get<double>("separator_width");
    this->collector_width = database->get<double>("collector_width");
    this->sandwich_height = database->get<double>("sandwich_height");
    this->tab_height      = database->get<double>("tab_height"     ); 

    std::string const mesh_file = database->get<std::string>("mesh_file");
    this->read_mesh(mesh_file);

    this->anode_tab_bbox = 
        std::make_pair(
            dealii::Point<dim>(-this->collector_width, this->sandwich_height), 
            dealii::Point<dim>(0.0, this->sandwich_height+this->tab_height)
        );
    this->anode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(-this->collector_width, 0.0), 
            dealii::Point<dim>(0.0, this->sandwich_height)
        );
    this->anode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(0.0, 0.0), 
            dealii::Point<dim>(this->electrode_width, this->sandwich_height)
        );
    this->separator_bbox =
        std::make_pair(
            dealii::Point<dim>(this->electrode_width, 0.0), 
            dealii::Point<dim>(this->electrode_width+this->separator_width, this->sandwich_height)
        );
    this->cathode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(this->electrode_width+this->separator_width, 0.0), 
            dealii::Point<dim>(this->electrode_width+this->separator_width+this->electrode_width, this->sandwich_height)
        );
    this->cathode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(this->electrode_width+this->separator_width+this->electrode_width, 0.0), 
            dealii::Point<dim>(this->electrode_width+this->separator_width+this->electrode_width+this->collector_width, this->sandwich_height)
        );
    this->cathode_tab_bbox =
        std::make_pair(
            dealii::Point<dim>(this->electrode_width+this->separator_width+this->electrode_width, -this->tab_height), 
            dealii::Point<dim>(this->electrode_width+this->separator_width+this->electrode_width+this->collector_width, 0.0)
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

    typename dealii::Triangulation<dim>::active_cell_iterator cell     = this->triangulation.begin_active();
    typename dealii::Triangulation<dim>::active_cell_iterator end_cell = this->triangulation.end();
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
 
//        for (std::size_t face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
//            if (cell->face(face)->at_boundary()) {
//                if ((std::abs(cell->face(face)->center()[1] - y_max) < 1.0e-10)
//                    && (cell->material_id() == anode_collector_material_id)) {
//                    cell->face(face)->set_boundary_indicator(anode_boundary_id);
//                } else if ((std::abs(cell->face(face)->center()[1] - y_min) < 1.0e-10)
//                    && (cell->material_id() == cathode_collector_material_id)) {
//                    cell->face(face)->set_boundary_indicator(cathode_boundary_id);
//                } else if ((std::abs(cell->face(face)->center()[1] - y_top) < 1.0e-10)
//                    && ((cell->material_id() == cathode_electrode_material_id)
//                        || (cell->material_id() == anode_electrode_material_id)
//                        || (cell->material_id() == separator_material_id))
//                    ) {
//                    cell->face(face)->set_boundary_indicator(upper_boundary_id);
//                } else if ((std::abs(cell->face(face)->center()[1] - y_bottom) < 1.0e-10)
//                    && ((cell->material_id() == cathode_electrode_material_id)
//                        || (cell->material_id() == anode_electrode_material_id)
//                        || (cell->material_id() == separator_material_id))
//                    ) {
//                    cell->face(face)->set_boundary_indicator(lower_boundary_id);
//                } else {
//                    cell->face(face)->set_boundary_indicator(other_boundary_id);
//                } // end if
//            } // end if face at boundary
//        } // end for face
    } // end for cell
    this->write_mesh("output_test_geometry_0.vtk");
}

template <int dim>
void
Geometry<dim>::
read_mesh(std::string const & mesh_file)
{
    dealii::GridIn<dim> mesh_reader;
    mesh_reader.attach_triangulation(this->triangulation);
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
void
Geometry<dim>::
write_mesh(std::string const & mesh_file) const
{
    dealii::GridOut mesh_writer;
    std::fstream fout;
    fout.open(mesh_file.c_str(), std::fstream::out);
    std::string const file_extension = mesh_file.substr(mesh_file.find_last_of(".")+1);
    if (file_extension.compare("vtk") == 0) {
        mesh_writer.write_vtk(this->triangulation, fout);
    } else {
        throw std::runtime_error("Bad output format ."+file_extension+" in mesh file "+mesh_file);
    }
    fout.close();
}

template <int dim>
void
Geometry<dim>::
move_mesh(std::shared_ptr<boost::property_tree::ptree const> const & database)
{
    this->electrode_width = database->get<double>("electrode_width");
    this->separator_width = database->get<double>("separator_width");
    this->collector_width = database->get<double>("collector_width");
    this->sandwich_height = database->get<double>("sandwich_height");
    this->tab_height      = database->get<double>("tab_height"     ); 
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_anode_tab_bbox = 
        std::make_pair(
            dealii::Point<dim>(0.0, this->sandwich_height), 
            dealii::Point<dim>(this->collector_width, this->sandwich_height+this->tab_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_anode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(0.0, 0.0), 
            dealii::Point<dim>(this->collector_width, this->sandwich_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_anode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(this->collector_width, 0.0), 
            dealii::Point<dim>(this->collector_width+this->electrode_width, this->sandwich_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_separator_bbox =
        std::make_pair(
            dealii::Point<dim>(this->collector_width+this->electrode_width, 0.0), 
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width, this->sandwich_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_cathode_electrode_bbox =
        std::make_pair(
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width, 0.0), 
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width+this->electrode_width, this->sandwich_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_cathode_collector_bbox =
        std::make_pair(
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width+this->electrode_width, 0.0), 
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width+this->electrode_width+this->collector_width, this->sandwich_height)
        );
    std::pair<dealii::Point<dim>,dealii::Point<dim> > new_cathode_tab_bbox =
        std::make_pair(
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width+this->electrode_width, -this->tab_height), 
            dealii::Point<dim>(this->collector_width+this->electrode_width+this->separator_width+this->electrode_width+this->collector_width, 0.0)
        );
    
    auto transform_point =
        [](dealii::Point<dim> const & p_0, dealii::Point<dim> const & p_1, dealii::Point<dim> const & p,
           dealii::Point<dim> const & q_0, dealii::Point<dim> const & q_1, dealii::Point<dim> &       q)
        {
            for (int d = 0; d < dim; ++d) {
                q[d] = (q_1[d] - q_0[d]) / (p_1[d] - p_0[d]) * (p[d] - p_0[d]) + q_0[d];
            } // end for d
        };
    std::vector<bool> vertex_visited(this->triangulation.n_vertices(), false);
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
            
    typename dealii::Triangulation<dim>::active_cell_iterator cell     = this->triangulation.begin_active();
    typename dealii::Triangulation<dim>::active_cell_iterator end_cell = this->triangulation.end();
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
    this->write_mesh("output_test_geometry_1.vtk");
}


} // end namespace cap

int main(int argc, char *argv[])
{
    std::shared_ptr<boost::property_tree::ptree> params(new boost::property_tree::ptree);
    params->put("electrode_width", 50.0e-6);
    params->put("separator_width", 25.0e-6);
    params->put("collector_width",  5.0e-6);
    params->put("sandwich_height", 25.0e-6);
    params->put("tab_height",       5.0e-6);
    params->put("separator_material_id",         3);
    params->put("anode_electrode_material_id",   4);
    params->put("anode_collector_material_id",   5);
    params->put("cathode_electrode_material_id", 6);
    params->put("cathode_collector_material_id", 7);
    params->put("mesh_file", "mesh_2d.ucd");

    cap::Geometry<2> geo(params);

    dealii::Triangulation<2> const & tria = geo.get_triangulation();
    std::cout<<"cells="<<tria.n_active_cells()<<"  "
        <<"faces="<<tria.n_active_faces()<<"  "
        <<"vertices="<<tria.n_used_vertices()<<"\n";

    params->put("electrode_width", 150.0e-6);
    params->put("separator_width",  25.0e-6);
    params->put("collector_width",  50.0e-6);
    params->put("sandwich_height",  75.0e-6);
    params->put("tab_height",        2.5e-6);
    geo.move_mesh(params);

    return 0;
}
