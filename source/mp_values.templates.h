#include <cache/mp_values.h>

#include <stdexcept>

template <int dim, int spacedim>
void
MPValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           std::vector<double> &        values) const

{
    throw std::runtime_error("Should have been overloaded..."); 
}

template <int dim, int spacedim>
void
MPValues<dim, spacedim>::
get_values(std::string const &                         key,
           active_cell_iterator const &                cell,
           std::vector<dealii::Tensor<1, spacedim> > & values) const

{
    throw std::runtime_error("Should have been overloaded..."); 
}

//////////////////////// SUPERCAPACITORS ////////////////////////////
template <int dim, int spacedim>
SuperCapacitorMPValues<dim, spacedim>::
SuperCapacitorMPValues(MPValuesParameters<dim, spacedim> const & parameters)
    : MPValues<dim, spacedim>(parameters)
{
    SuperCapacitorMPValuesParameters<dim, spacedim> const * super_capacitor_parameters = 
        dynamic_cast<SuperCapacitorMPValuesParameters<dim, spacedim> const *>(&parameters);

    this->anode_collector_material_id   = super_capacitor_parameters->anode_collector_material_id;
    this->cathode_collector_material_id = super_capacitor_parameters->cathode_collector_material_id;
    this->anode_material_id             = super_capacitor_parameters->anode_material_id;
    this->cathode_material_id           = super_capacitor_parameters->cathode_material_id;
    this->separator_material_id         = super_capacitor_parameters->separator_material_id;

    this->separator_thermal_conductivity = super_capacitor_parameters->separator_thermal_conductivity; 
    this->electrode_thermal_conductivity = super_capacitor_parameters->electrode_thermal_conductivity; 
    this->collector_thermal_conductivity = super_capacitor_parameters->collector_thermal_conductivity; 
    this->separator_density              = super_capacitor_parameters->separator_density;              
    this->electrode_density              = super_capacitor_parameters->electrode_density;              
    this->collector_density              = super_capacitor_parameters->collector_density;              
    this->separator_heat_capacity        = super_capacitor_parameters->separator_heat_capacity;        
    this->electrode_heat_capacity        = super_capacitor_parameters->electrode_heat_capacity;        
    this->collector_heat_capacity        = super_capacitor_parameters->collector_heat_capacity;        

    this->specific_capacitance           = super_capacitor_parameters->specific_capacitance;          
    this->electrode_void_volume_fraction = super_capacitor_parameters->electrode_void_volume_fraction;
    this->separator_void_volume_fraction = super_capacitor_parameters->separator_void_volume_fraction;
    this->electrolyte_conductivity       = super_capacitor_parameters->electrolyte_conductivity;      
    this->solid_phase_conductivity       = super_capacitor_parameters->solid_phase_conductivity;      

}

template <int dim, int spacedim>
void 
SuperCapacitorMPValues<dim, spacedim>::
get_values(std::string const &          key,
           active_cell_iterator const & cell,
           std::vector<double> &        values) const
{
    if (key.compare("thermal_conductivity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), this->separator_thermal_conductivity);
        } else if ((cell->material_id() == this->anode_material_id) 
            || (cell->material_id() == this->cathode_material_id)) {
            std::fill(values.begin(), values.end(), this->electrode_thermal_conductivity);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), this->collector_thermal_conductivity);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("density_times_heat_capacity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), this->separator_density*this->separator_heat_capacity);
        } else if ((cell->material_id() == this->anode_material_id) 
            || (cell->material_id() == this->cathode_material_id)) {
            std::fill(values.begin(), values.end(), this->electrode_density*this->electrode_heat_capacity);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), this->collector_density*this->collector_heat_capacity);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("specific_capacitance") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if ((cell->material_id() == this->anode_material_id) 
            || (cell->material_id() == this->cathode_material_id)) {
            std::fill(values.begin(), values.end(), this->specific_capacitance);
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("solid_electrical_conductivity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), 0.0);
        } else if ((cell->material_id() == this->anode_material_id) 
            || (cell->material_id() == this->cathode_material_id)) {
            std::fill(values.begin(), values.end(), this->solid_phase_conductivity*std::pow((1.0-this->electrode_void_volume_fraction), 1.5));
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), this->solid_phase_conductivity);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else if (key.compare("liquid_electrical_conductivity") == 0) {
        if (cell->material_id() == this->separator_material_id) {
            std::fill(values.begin(), values.end(), this->electrolyte_conductivity*std::pow(this->separator_void_volume_fraction, 1.5));
        } else if ((cell->material_id() == this->anode_material_id) 
            || (cell->material_id() == this->cathode_material_id)) {
            std::fill(values.begin(), values.end(), this->electrolyte_conductivity*std::pow(this->electrode_void_volume_fraction, 1.5));
        } else if ((cell->material_id() == this->anode_collector_material_id) 
            || (cell->material_id() == this->cathode_collector_material_id)) {
            std::fill(values.begin(), values.end(), 0.0);
        } else {
            throw std::runtime_error("Invalid material id");
        } // end if material id
    } else {
        throw std::runtime_error("Invalid key");
    } // end if key
}
