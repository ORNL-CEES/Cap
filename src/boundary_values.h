#ifndef CAP_BOUNDARY_VALUES_H
#define CAP_BOUNDARY_VALUES_H

//#include <deal.II/base/types.h>
//#include <deal.II/base/tensor.h>
//#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <boost/property_tree/ptree.hpp>

namespace cap {

//////////////////////// BOUNDARY VALUES PARAMETERS ////////////////////////////
template <int dim, int spacedim=dim>
class BoundaryValuesParameters {
public:
    BoundaryValuesParameters(std::shared_ptr<boost::property_tree::ptree> d) 
        : database(d)
    { }
    virtual ~BoundaryValuesParameters() { }
// keep public for now
    std::shared_ptr<boost::property_tree::ptree> database;

};

//////////////////////// BOUNDARY VALUES ////////////////////////////
template <int dim, int spacedim=dim>
class BoundaryValues {
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    BoundaryValues(BoundaryValuesParameters<dim, spacedim> const & ) { }
    virtual ~BoundaryValues() { }
    virtual void get_values(std::string const &          key,
                            active_cell_iterator const & cell,
                            unsigned int const           face,
                            std::vector<double> &        values) const;
    virtual void get_values(std::string const &                         key,
                            active_cell_iterator const &                cell,
                            unsigned int const                          face,
                            std::vector<dealii::Tensor<1, spacedim> > & values) const;
};


//////////////////////// SUPERCAPACITOR BOUNDARY VALUES ////////////////////////////
template <int dim, int spacedim=dim>
class SuperCapacitorBoundaryValuesParameters : public BoundaryValuesParameters<dim, spacedim> {
public:
    SuperCapacitorBoundaryValuesParameters(std::shared_ptr<boost::property_tree::ptree> d)
        : BoundaryValuesParameters<dim, spacedim>(d)
    { }
};

template <int dim, int spacedim=dim>
class SuperCapacitorBoundaryValues : public BoundaryValues<dim, spacedim> { 
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    SuperCapacitorBoundaryValues(BoundaryValuesParameters<dim, spacedim> const & parameters);
    void get_values(std::string const &          key,
                    active_cell_iterator const & cell,
                    unsigned int const           face,
                    std::vector<double> &        values) const override;
protected:
    dealii::types::material_id separator_material_id;
    dealii::types::material_id anode_electrode_material_id;
    dealii::types::material_id anode_collector_material_id;
    dealii::types::material_id cathode_electrode_material_id;
    dealii::types::material_id cathode_collector_material_id;

    dealii::types::boundary_id anode_boundary_id;
    dealii::types::boundary_id cathode_boundary_id; 
    dealii::types::boundary_id upper_boundary_id;
    dealii::types::boundary_id lower_boundary_id;
    dealii::types::boundary_id other_boundary_id;

    double charge_potential; 
    double discharge_potential; 
    double charge_current_density; 
    double discharge_current_density;
    double initial_potential;

    double upper_ambient_temperature;
    double lower_ambient_temperature;
    double upper_heat_transfer_coefficient;
    double lower_heat_transfer_coefficient;
};

} // end namespace cap

#endif // CAP_BOUNDARY_VALUES_H
