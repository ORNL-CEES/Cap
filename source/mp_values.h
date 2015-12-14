#ifndef CAP_MP_VALUES_H
#define CAP_MP_VALUES_H

#include <cap/geometry.h>
//#include <deal.II/base/types.h>
//#include <deal.II/base/tensor.h>
//#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <boost/property_tree/ptree.hpp>

namespace cap {

//////////////////////// MP VALUES PARAMETERS ////////////////////////////
template <int dim, int spacedim=dim>
class MPValuesParameters {
public:
    MPValuesParameters(std::shared_ptr<boost::property_tree::ptree const> d)
        : database(d)
    { }
    virtual ~MPValuesParameters() = default;
// keep public for now
    std::shared_ptr<boost::property_tree::ptree const> database;
    std::shared_ptr<Geometry<dim> const> geometry;

};

//////////////////////// MP VALUES ////////////////////////////
template <int dim, int spacedim=dim>
class MPValues {
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    MPValues(MPValuesParameters<dim, spacedim> const & params);
    virtual ~MPValues() = default;
    virtual void get_values(std::string const &          key,
                            active_cell_iterator const & cell,
                            std::vector<double> &        values) const;
    virtual void get_values(std::string const &                         key,
                            active_cell_iterator const &                cell,
                            std::vector<dealii::Tensor<1, spacedim> > & values) const;
protected:
    std::unordered_map<dealii::types::material_id, std::shared_ptr<MPValues<dim>>> materials;
};


//////////////////////// SUPERCAPACITOR MP VALUES ////////////////////////////
template <int dim, int spacedim=dim>
class SuperCapacitorMPValuesParameters : public MPValuesParameters<dim, spacedim> {
public:
    SuperCapacitorMPValuesParameters(std::shared_ptr<boost::property_tree::ptree const> d)
        : MPValuesParameters<dim, spacedim>(d)
    { }
};

template <int dim, int spacedim=dim>
class SuperCapacitorMPValues : public MPValues<dim, spacedim> { 
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    SuperCapacitorMPValues(MPValuesParameters<dim, spacedim> const & parameters);
    void get_values(std::string const &          key,
                    active_cell_iterator const & cell,
                    std::vector<double> &        values) const override;
protected:
    dealii::types::material_id separator_material_id;
    dealii::types::material_id anode_electrode_material_id;
    dealii::types::material_id anode_collector_material_id;
    dealii::types::material_id cathode_electrode_material_id;
    dealii::types::material_id cathode_collector_material_id;

    double separator_thermal_conductivity;
    double electrode_thermal_conductivity;
    double collector_thermal_conductivity;
    double separator_density;
    double electrode_density;
    double collector_density;
    double electrolyte_mass_density;
    double separator_heat_capacity;
    double electrode_heat_capacity;
    double collector_heat_capacity;

    double differential_capacitance;
    double electrode_void_volume_fraction;
    double separator_void_volume_fraction;
    double electrolyte_conductivity;
    double solid_phase_conductivity;
    double collector_electrical_resistivity;
    double separator_tortuosity_factor;
    double electrode_tortuosity_factor;
    double bruggemans_coefficient;
    double pores_characteristic_dimension;
    double pores_geometry_factor;
    double specific_surface_area_per_unit_volume;

    double anodic_charge_transfer_coefficient;
    double cathodic_charge_transfer_coefficient;
    double faraday_constant;
    double gas_constant;
    double exchange_current_density;
    double temperature;
};



//////////////////////// NEW STUFF ////////////////////////////
template <int dim, int spacedim=dim>
class NewStuffMPValues : public MPValues<dim, spacedim> {
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    NewStuffMPValues(MPValuesParameters<dim, spacedim> const & parameters);
    void get_values(std::string const &          key,
                    active_cell_iterator const & cell,
                    std::vector<double> &        values) const override;
protected:
    std::unordered_map<std::string, std::function<void(active_cell_iterator const &, std::vector<double> &)>> properties;
};



template <int dim, int spacedim=dim>
class PorousElectrodeMPValues : public NewStuffMPValues<dim, spacedim> {
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    PorousElectrodeMPValues(MPValuesParameters<dim, spacedim> const & parameters);
};



template <int dim, int spacedim=dim>
class MetalFoilMPValues : public NewStuffMPValues<dim, spacedim> {
public:
    typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator active_cell_iterator;
    MetalFoilMPValues(MPValuesParameters<dim, spacedim> const & parameters);
};



template <int dim, int spacedim=dim>
std::shared_ptr<MPValues<dim>>
buildMaterial(std::string const & material_name, std::shared_ptr<boost::property_tree::ptree const> database)
{
    std::shared_ptr<boost::property_tree::ptree> material_database =
        std::make_shared<boost::property_tree::ptree>(database->get_child(material_name));
    std::string const type = material_database->get<std::string>("type");
    if (type.compare("deprecated_stuff") == 0) {
        return std::make_shared<SuperCapacitorMPValues<dim>>(MPValuesParameters<dim>(database));
    } else if (type.compare("porous_electrode") == 0) {
        std::shared_ptr<boost::property_tree::ptree> dummy_database =
            std::make_shared<boost::property_tree::ptree>(*database);
        dummy_database->put("ugly_hack", material_name);
        return std::make_shared<PorousElectrodeMPValues<dim>>(MPValuesParameters<dim>(dummy_database));
    } else if (type.compare("permeable_membrane") == 0) {
        std::shared_ptr<boost::property_tree::ptree> dummy_database =
            std::make_shared<boost::property_tree::ptree>(*database);
        dummy_database->put("ugly_hack", material_name);
        std::string const matrix_phase = dummy_database->get<std::string>(material_name+"."+"matrix_phase");
        dummy_database->put(matrix_phase+"."+"differential_capacitance", 0.0);
        dummy_database->put(matrix_phase+"."+"exchange_current_density", 0.0);
        dummy_database->put(matrix_phase+"."+"electrical_conductivity" , 0.0);
        return std::make_shared<PorousElectrodeMPValues<dim>>(MPValuesParameters<dim>(dummy_database));
    } else if (type.compare("current_collector") == 0) {
        std::shared_ptr<boost::property_tree::ptree> dummy_database =
            std::make_shared<boost::property_tree::ptree>(*database);
        dummy_database->put("ugly_hack", material_name);
        return std::make_shared<MetalFoilMPValues<dim>>(MPValuesParameters<dim>(dummy_database));
    } else {
        throw std::runtime_error("Invalid material type "+type);
    }
}

} // end namespace cap

#endif // CAP_MP_VALUES_H
