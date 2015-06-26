#define BOOST_TEST_MODULE TestGeometry
#define BOOST_TEST_MAIN
#include <cap/utils.h>
#include <cap/geometry.h>
#include <deal.II/grid/grid_out.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>

template <int dim>
void
write_mesh(std::string const & mesh_file, std::shared_ptr<dealii::Triangulation<dim> const> triangulation)
{
    dealii::GridOut mesh_writer;
    std::fstream fout;
    fout.open(mesh_file.c_str(), std::fstream::out);
    std::string const file_extension = mesh_file.substr(mesh_file.find_last_of(".")+1);
    if (file_extension.compare("vtk") == 0) {
        mesh_writer.write_vtk(*triangulation, fout);
    } else {
        throw std::runtime_error("Bad output format ."+file_extension+" in mesh file "+mesh_file);
    }
    fout.close();
}

BOOST_AUTO_TEST_CASE( test_reset_geometry )
{
    std::shared_ptr<boost::property_tree::ptree> params(new boost::property_tree::ptree);
    params->put("electrode_width"              ,  50.0e-6     );
    params->put("separator_width"              ,  25.0e-6     );
    params->put("collector_width"              ,   5.0e-6     );
    params->put("sandwich_height"              ,  25.0e-6     );
    params->put("tab_height"                   ,   5.0e-6     );
    params->put("separator_material_id"        ,   3          );
    params->put("anode_electrode_material_id"  ,   4          );
    params->put("anode_collector_material_id"  ,   5          );
    params->put("cathode_electrode_material_id",   6          );
    params->put("cathode_collector_material_id",   7          );
    params->put("mesh_file"                    , "mesh_2d.ucd");
    params->put("materials"                    ,   2          );
    params->put("material_0.name"              , "electrode"  );
    params->put("material_0.material_id"       , "4,6"        );
    params->put("material_1.name"              , "separator"  );
    params->put("material_1.material_id"       , "3"          );

    cap::SuperCapacitorGeometry<2> geo(params);
    write_mesh("output_test_geometry_0.vtk", geo.get_triangulation());

    dealii::Triangulation<2> const & tria = *geo.get_triangulation();
    std::cout<<"cells="<<tria.n_active_cells()<<"  "
        <<"faces="<<tria.n_active_faces()<<"  "
        <<"vertices="<<tria.n_used_vertices()<<"\n";

    params->put("electrode_width"              , 150.0e-6     );
    params->put("separator_width"              ,  25.0e-6     );
    params->put("collector_width"              ,  50.0e-6     );
    params->put("sandwich_height"              ,  75.0e-6     );
    params->put("tab_height"                   ,   2.5e-6     );

    geo.reset(params);
    write_mesh("output_test_geometry_1.vtk", geo.get_triangulation());
}
