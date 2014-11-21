#include <cache/boundary_values.h>

#include <stdexcept>

template <int dim, int spacedim>
void
BoundaryValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           unsigned int const           face,
           std::vector<double> &        values) const

{
    throw std::runtime_error("Should have been overloaded..."); 
}

template <int dim, int spacedim>
void
BoundaryValues<dim, spacedim>::
get_values(std::string const &                         key,
           active_cell_iterator const &                cell,
           unsigned int const                          face,
           std::vector<dealii::Tensor<1, spacedim> > & values) const

{
    throw std::runtime_error("Should have been overloaded..."); 
}

//////////////////////// SUPERCAPACITORS ////////////////////////////
template <int dim, int spacedim>
SuperCapacitorBoundaryValues<dim, spacedim>::
SuperCapacitorBoundaryValues(BoundaryValuesParameters<dim, spacedim> const & parameters)
    : BoundaryValues<dim, spacedim>(parameters)
{
    SuperCapacitorBoundaryValuesParameters<dim, spacedim> const * super_capacitor_parameters = 
        dynamic_cast<SuperCapacitorBoundaryValuesParameters<dim, spacedim> const *>(&parameters);

    this->anode_material_id             = super_capacitor_parameters->anode_material_id;
    this->cathode_material_id           = super_capacitor_parameters->cathode_material_id;
    this->separator_material_id         = super_capacitor_parameters->separator_material_id;
    this->anode_collector_material_id   = super_capacitor_parameters->anode_collector_material_id;
    this->cathode_collector_material_id = super_capacitor_parameters->cathode_collector_material_id;

    this->charge_potential          = super_capacitor_parameters->charge_potential;        
    this->discharge_potential       = super_capacitor_parameters->discharge_potential;     
    this->charge_current_density    = super_capacitor_parameters->charge_current_density;  
    this->discharge_current_density = super_capacitor_parameters->discharge_current_density;
    this->initial_potential         = super_capacitor_parameters->initial_potential;   

    this->upper_ambient_temperature       = super_capacitor_parameters->upper_ambient_temperature; 
    this->lower_ambient_temperature       = super_capacitor_parameters->lower_ambient_temperature; 
    this->upper_heat_transfer_coefficient = super_capacitor_parameters->upper_heat_transfer_coefficient; 
    this->lower_heat_transfer_coefficient = super_capacitor_parameters->lower_heat_transfer_coefficient; 

    this->cathode_boundary_id = super_capacitor_parameters->cathode_boundary_id;
    this->anode_boundary_id   = super_capacitor_parameters->anode_boundary_id;          
    this->upper_boundary_id   = super_capacitor_parameters->upper_boundary_id;
    this->lower_boundary_id   = super_capacitor_parameters->lower_boundary_id;      
    this->other_boundary_id   = super_capacitor_parameters->other_boundary_id;      

}

template <int dim, int spacedim>
void 
SuperCapacitorBoundaryValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           unsigned int const           face,
           std::vector<double> &        values) const
{
    if (key.compare("charge_current_density") == 0) {
        if (cell->face(face)->boundary_indicator() == this->anode_boundary_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if (cell->face(face)->boundary_indicator() == this->cathode_boundary_id) {
            std::fill(values.begin(), values.end(), this->charge_current_density);
        } else if ((cell->face(face)->boundary_indicator() == this->upper_boundary_id) 
            || (cell->face(face)->boundary_indicator() == this->lower_boundary_id)
            || (cell->face(face)->boundary_indicator() == this->other_boundary_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid boundary id");
        } // end if boundary id
    } else if (key.compare("discharge_current_density") == 0) {
        if (cell->face(face)->boundary_indicator() == this->anode_boundary_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if (cell->face(face)->boundary_indicator() == this->cathode_boundary_id) {
            std::fill(values.begin(), values.end(), this->discharge_current_density);
        } else if ((cell->face(face)->boundary_indicator() == this->upper_boundary_id) 
            || (cell->face(face)->boundary_indicator() == this->lower_boundary_id)
            || (cell->face(face)->boundary_indicator() == this->other_boundary_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid boundary id");
        } // end if boundary id
    } else if (key.compare("ambient_temperature") == 0) {
        if ((cell->face(face)->boundary_indicator() == this->anode_boundary_id)
            || (cell->face(face)->boundary_indicator() == this->cathode_boundary_id)
            || (cell->face(face)->boundary_indicator() == this->other_boundary_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if (cell->face(face)->boundary_indicator() == this->upper_boundary_id) {
            std::fill(values.begin(), values.end(), this->upper_ambient_temperature);
        } else if (cell->face(face)->boundary_indicator() == this->lower_boundary_id) {
            std::fill(values.begin(), values.end(), this->lower_ambient_temperature);
        } else {
            throw std::runtime_error("Invalid boundary id");
        } // end if boundary id
    } else if (key.compare("heat_transfer_coefficient") == 0) {
        if ((cell->face(face)->boundary_indicator() == this->anode_boundary_id)
            || (cell->face(face)->boundary_indicator() == this->cathode_boundary_id)
            || (cell->face(face)->boundary_indicator() == this->other_boundary_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if (cell->face(face)->boundary_indicator() == this->upper_boundary_id) {
            std::fill(values.begin(), values.end(), this->upper_heat_transfer_coefficient);
        } else if (cell->face(face)->boundary_indicator() == this->lower_boundary_id) {
            std::fill(values.begin(), values.end(), this->lower_heat_transfer_coefficient);
        } else {
            throw std::runtime_error("Invalid boundary id");
        } // end if boundary id
    } else {
        throw std::runtime_error("Invalid key");
    } // end if key
}
