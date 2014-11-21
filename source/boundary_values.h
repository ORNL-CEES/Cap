#ifndef CACHE_BOUNDARY_VALUES_H
#define CACHE_BOUNDARY_VALUES_H

//#include <deal.II/base/types.h>
//#include <deal.II/base/tensor.h>
//#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

//////////////////////// BOUNDARY VALUES PARAMETERS ////////////////////////////
template <int dim, int spacedim=dim>
class BoundaryValuesParameters {
public:
    virtual ~BoundaryValuesParameters() { }

};

//////////////////////// BOUNDARY VALUES ////////////////////////////
template <int dim, int spacedim=dim>
class BoundaryValues {
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    BoundaryValues(BoundaryValuesParameters<dim, spacedim> const & parameters) { }
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
    dealii::types::material_id anode_collector_material_id;
    dealii::types::material_id cathode_collector_material_id;
    dealii::types::material_id separator_material_id;
    dealii::types::material_id anode_material_id;
    dealii::types::material_id cathode_material_id;

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

template <int dim, int spacedim=dim>
class SuperCapacitorBoundaryValues : public BoundaryValues<dim, spacedim> { 
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    SuperCapacitorBoundaryValues(BoundaryValuesParameters<dim, spacedim> const & parameters);
    void get_values(std::string const &          key,
                    active_cell_iterator const & cell,
                    unsigned int const           face,
                    std::vector<double> &        values) const;
protected:
    dealii::types::material_id anode_collector_material_id;
    dealii::types::material_id cathode_collector_material_id;
    dealii::types::material_id separator_material_id;
    dealii::types::material_id anode_material_id;
    dealii::types::material_id cathode_material_id;

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

#endif // CACHE_BOUNDARY_VALUES_H
