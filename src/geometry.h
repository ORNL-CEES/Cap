#ifndef CAP_GEOMETRY_H
#define CAP_GEOMETRY_H

#include <deal.II/base/types.h>
#include <deal.II/grid/tria.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace cap {

template <int dim>
class Geometry
{
public:
    Geometry(std::shared_ptr<boost::property_tree::ptree const> const & database);
    virtual ~Geometry() = default;
    inline std::shared_ptr<dealii::Triangulation<dim> const> get_triangulation() const
        { return this->triangulation; }
    virtual void reset(std::shared_ptr<boost::property_tree::ptree const> const & database) = 0;
protected:
    std::shared_ptr<dealii::Triangulation<dim> > triangulation;
};



template <int dim>
class SuperCapacitorGeometry : public Geometry<dim>
{
public:
    SuperCapacitorGeometry(std::shared_ptr<boost::property_tree::ptree const> const & database);
    void reset(std::shared_ptr<boost::property_tree::ptree const> const & database) override;

private:
    dealii::types::material_id separator_material_id        ;
    dealii::types::material_id anode_electrode_material_id  ;
    dealii::types::material_id anode_collector_material_id  ;
    dealii::types::material_id cathode_electrode_material_id;
    dealii::types::material_id cathode_collector_material_id;

    std::pair<dealii::Point<dim>,dealii::Point<dim> > anode_tab_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > anode_collector_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > anode_electrode_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > separator_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > cathode_electrode_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > cathode_collector_bbox;
    std::pair<dealii::Point<dim>,dealii::Point<dim> > cathode_tab_bbox;
    static int const spacedim = dim;
};

} // end namespace cap

#endif // CAP_GEOMETRY_H
