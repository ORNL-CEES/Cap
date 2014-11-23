#include <cache/electrochemical_operator.h>
#include <cache/dof_extractor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>

namespace cache {

template <int dim>
ElectrochemicalOperator<dim>::
ElectrochemicalOperator(std::shared_ptr<OperatorParameters<dim> const> parameters)
  : Operator<dim>(parameters) 
{
    std::shared_ptr<boost::property_tree::ptree const> database = parameters->database;

    this->solid_potential_component  = database->get<unsigned int>("solid_potential_component" );
    this->liquid_potential_component = database->get<unsigned int>("liquid_potential_component");
    this->temperature_component      = database->get<unsigned int>("temperature_component"     );
                                              
    this->charge_potential           = database->get<double>("boundary_values.charge_potential"         );
    this->discharge_potential        = database->get<double>("boundary_values.discharge_potential"      );
    this->charge_current_density     = database->get<double>("boundary_values.charge_current_density"   );
    this->discharge_current_density  = database->get<double>("boundary_values.discharge_current_density");
    this->initial_potential          = database->get<double>("boundary_values.initial_potential"        );
                                              
    this->alpha                      = database->get<double>("material_properties.alpha");
                                              
    this->anode_boundary_id          = database->get<dealii::types::boundary_id>("boundary_values.anode_boundary_id"  );
    this->cathode_boundary_id        = database->get<dealii::types::boundary_id>("boundary_values.cathode_boundary_id");

    std::shared_ptr<ElectrochemicalOperatorParameters<dim> const> electrochemical_parameters = 
        std::dynamic_pointer_cast<ElectrochemicalOperatorParameters<dim> const>(parameters);
}

template <int dim>
void
ElectrochemicalOperator<dim>::
reset(std::shared_ptr<OperatorParameters<dim> const> parameters)
{
std::cout<<"### reset electrical ###\n";
    std::shared_ptr<ElectrochemicalOperatorParameters<dim> const> electrochemical_parameters = 
        std::dynamic_pointer_cast<ElectrochemicalOperatorParameters<dim> const>(parameters);

    this->capacitor_state = electrochemical_parameters->capacitor_state;

    this->stiffness_matrix = 0.0;
    this->mass_matrix = 0.0;
    this->load_vector = 0.0;
    this->boundary_values.clear();

    this->compute_electrical_operator_contribution();

    if ((this->capacitor_state == GalvanostaticCharge)
        || (this->capacitor_state == GalvanostaticDischarge)
        || (this->capacitor_state == Relaxation)) {
        this->compute_neumann_boundary_contribution();
    } else if ((this->capacitor_state == PotentiostaticCharge)
        || (this->capacitor_state == PotentiostaticDischarge)
        || (this->capacitor_state == Initialize)) {
        this->compute_dirichlet_boundary_values();
    } else {
        throw std::runtime_error("State of the capacitor undetermined");
    } // end if


std::cout<<std::setprecision(15)
    <<"stiffness="<<this->stiffness_matrix.l1_norm()<<"  "
    <<"mass="<<this->mass_matrix.l1_norm()<<"  "
    <<"load="<<this->load_vector.l2_norm()<<"\n";
}

template <int dim>
void 
ElectrochemicalOperator<dim>::
compute_dirichlet_boundary_values() {
double potential;
if (this->capacitor_state == PotentiostaticCharge) {
    potential = this->charge_potential;
} else if (this->capacitor_state == PotentiostaticDischarge) {
    potential = this->discharge_potential;
} else if (this->capacitor_state == Initialize) {
    potential = this->initial_potential;
} else {
    throw std::runtime_error("What are you doing here?");
}
    
    unsigned int const n_components = dealii::DoFTools::n_components(this->dof_handler);
    std::vector<bool> mask(n_components, false);
    mask[this->solid_potential_component] = true;
    dealii::ComponentMask component_mask(mask);
    typename dealii::FunctionMap<dim>::type dirichlet_boundary_condition;
    std::vector<dealii::std_cxx11::shared_ptr<dealii::Function<dim> > > dirichlet_functions;
    dirichlet_functions.push_back(dealii::std_cxx11::shared_ptr<dealii::Function<dim> >(new dealii::ConstantFunction<dim>(0.0, n_components)));
    dirichlet_functions.push_back(dealii::std_cxx11::shared_ptr<dealii::Function<dim> >(new dealii::ConstantFunction<dim>(potential, n_components)));
    dirichlet_boundary_condition[this->anode_boundary_id] = dirichlet_functions[0].get();
    dirichlet_boundary_condition[this->cathode_boundary_id] = dirichlet_functions[1].get();
    dealii::VectorTools::interpolate_boundary_values(this->dof_handler, dirichlet_boundary_condition, this->boundary_values, component_mask);
}

template <int dim>
void 
ElectrochemicalOperator<dim>::
compute_neumann_boundary_contribution() {
{
    unsigned int const n_components = dealii::DoFTools::n_components(this->dof_handler);
    std::vector<bool> mask(n_components, false);
    mask[this->solid_potential_component] = true;
    dealii::ComponentMask component_mask(mask);
    dealii::VectorTools::interpolate_boundary_values(this->dof_handler, this->anode_boundary_id, dealii::ConstantFunction<dim>(0.0, n_components), this->boundary_values, component_mask);
}

std::string current_density;
if (this->capacitor_state == GalvanostaticCharge) {
    current_density = "charge_current_density";
} else if (this->capacitor_state == GalvanostaticDischarge) {
    current_density = "discharge_current_density";
} else if (this->capacitor_state == Relaxation) {
    return;
} else {
    throw std::runtime_error("Did you get lost?");
}

    dealii::FEValuesExtractors::Scalar const solid_potential(this->solid_potential_component);
    dealii::FiniteElement<dim> const & fe = this->dof_handler.get_fe();
    dealii::QGauss<dim-1> face_quadrature_rule(fe.degree+1); // TODO: maybe use fe extractor
    dealii::FEFaceValues<dim> fe_face_values(fe, face_quadrature_rule,
        dealii::update_values | dealii::update_JxW_values);
    unsigned int const dofs_per_cell = fe.dofs_per_cell;
    unsigned int const n_face_q_points = face_quadrature_rule.size();
    dealii::Vector<double> cell_load_vector(dofs_per_cell);
    std::vector<double> current_density_values(n_face_q_points);
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
unsigned int const n_components = dealii::DoFTools::n_components(this->dof_handler);
dealii::ComponentMask mask(n_components, false);
mask.set(this->solid_potential_component, true);
mask.set(this->liquid_potential_component, true);
DoFExtractor dof_extractor(mask, mask, dofs_per_cell);
    typename dealii::DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        end_cell = this->dof_handler.end();
    for ( ; cell != end_cell; ++cell) {
        cell_load_vector = 0.0;
        if (cell->at_boundary()) {
            for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {
                if (cell->face(face)->at_boundary()) {
                    fe_face_values.reinit(cell, face);
                    (this->b_values)->get_values(current_density, cell, face, current_density_values);
                    for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) {
                        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                            cell_load_vector(i) += (
                                  current_density_values[q_point] *
                                  fe_face_values[solid_potential].value(i, q_point)
                                ) * fe_face_values.JxW(q_point);
                        } // end for i
                    } // end for face quadrature point
                } // end if face at boundary
            } // end for face
        } // end if cell at boundary
        cell->get_dof_indices(local_dof_indices);
std::vector<dealii::types::global_dof_index> tmp_indices = dof_extractor.extract_row_indices(local_dof_indices);        
dealii::Vector<double> tmp_load_vector = dof_extractor.extract_vector(cell_load_vector);
        this->constraint_matrix.distribute_local_to_global(tmp_load_vector, tmp_indices, this->load_vector);
//        this->constraint_matrix.distribute_local_to_global(cell_load_vector, local_dof_indices, this->load_vector);
    } // end for cell

}
template <int dim>
void 
ElectrochemicalOperator<dim>::
compute_electrical_operator_contribution() 
{
    dealii::FiniteElement<dim> const & fe = this->dof_handler.get_fe();
    dealii::FEValuesExtractors::Scalar const solid_potential(this->solid_potential_component);
    dealii::FEValuesExtractors::Scalar const liquid_potential(this->liquid_potential_component);
    dealii::QGauss<dim> quadrature_rule(fe.degree+1); // TODO: map to component...
    dealii::FEValues<dim> fe_values(fe, quadrature_rule, 
        dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);
    unsigned int const dofs_per_cell = fe.dofs_per_cell;
    unsigned int const n_q_points = quadrature_rule.size();
    dealii::FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
    dealii::FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<double> solid_phase_diffusion_coefficient_values(n_q_points);
    std::vector<double> liquid_phase_diffusion_coefficient_values(n_q_points);
    std::vector<double> specific_capacitance_values(n_q_points);
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
unsigned int const n_components = dealii::DoFTools::n_components(this->dof_handler);
dealii::ComponentMask mask(n_components, false);
mask.set(this->solid_potential_component, true);
mask.set(this->liquid_potential_component, true);
DoFExtractor dof_extractor(mask, mask, dofs_per_cell);
    typename dealii::DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        end_cell = this->dof_handler.end();
    for ( ; cell != end_cell; ++cell) {
        cell_stiffness_matrix = 0.0;
        cell_mass_matrix = 0.0;
        fe_values.reinit(cell);
        (this->mp_values)->get_values("specific_capacitance", cell, specific_capacitance_values);
        (this->mp_values)->get_values("solid_electrical_conductivity", cell, solid_phase_diffusion_coefficient_values);
        (this->mp_values)->get_values("liquid_electrical_conductivity", cell, liquid_phase_diffusion_coefficient_values);
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    cell_stiffness_matrix(i, j) += (
                        + solid_phase_diffusion_coefficient_values[q_point] *
                        ( fe_values[solid_potential].gradient(i, q_point) * 
                          fe_values[solid_potential].gradient(j, q_point) )
                        + liquid_phase_diffusion_coefficient_values[q_point] *
                        ( fe_values[liquid_potential].gradient(i, q_point) * 
                          fe_values[liquid_potential].gradient(j, q_point) )
                        ) * fe_values.JxW(q_point);
                    cell_mass_matrix(i, j) += (
                        + specific_capacitance_values[q_point] *
                        ( fe_values[solid_potential].value(i, q_point) * 
                          fe_values[solid_potential].value(j, q_point) )
                        - specific_capacitance_values[q_point] *
                        ( fe_values[solid_potential].value(i, q_point) * 
                          fe_values[liquid_potential].value(j, q_point) ) 
                        - specific_capacitance_values[q_point] *
                        ( fe_values[liquid_potential].value(i, q_point) * 
                          fe_values[solid_potential].value(j, q_point) )
                        + specific_capacitance_values[q_point] *
                        ( fe_values[liquid_potential].value(i, q_point) * 
                          fe_values[liquid_potential].value(j, q_point) ) 
                        ) * fe_values.JxW(q_point);
                } // end for j
            } // end for i
        } // end for quadrature point
        cell->get_dof_indices(local_dof_indices);
std::vector<dealii::types::global_dof_index> tmp_indices = dof_extractor.extract_row_indices(local_dof_indices);        
dealii::FullMatrix<double> tmp_mass_matrix = dof_extractor.extract_matrix(cell_mass_matrix);
dealii::FullMatrix<double> tmp_stiffness_matrix = dof_extractor.extract_matrix(cell_stiffness_matrix);
        this->constraint_matrix.distribute_local_to_global(tmp_stiffness_matrix, tmp_indices, this->stiffness_matrix);
        this->constraint_matrix.distribute_local_to_global(tmp_mass_matrix, tmp_indices, this->mass_matrix);

//        this->constraint_matrix.distribute_local_to_global(cell_stiffness_matrix, local_dof_indices, this->stiffness_matrix);
//        this->constraint_matrix.distribute_local_to_global(cell_mass_matrix, local_dof_indices, this->mass_matrix);
    } // end for cell
}
template <int dim>
void 
ElectrochemicalOperator<dim>::
compute_heat_source(dealii::BlockVector<double> const & potential_solution_vector,
                    dealii::Vector<double> &            thermal_load_vector) const
{
    double coeff;
    if (this->capacitor_state == GalvanostaticCharge) {
        coeff = this->alpha;
    } else if (this->capacitor_state == GalvanostaticDischarge) {
        coeff = - this->alpha;
    } else if (this->capacitor_state == Relaxation) {
        coeff = 0.0;
    } else {
        throw std::runtime_error("compute_heat_source irreversible component not implemented yet");
    } // end if
    dealii::FiniteElement<dim> const & fe = this->dof_handler.get_fe();
    dealii::FEValuesExtractors::Scalar const solid_potential(this->solid_potential_component);
    dealii::FEValuesExtractors::Scalar const liquid_potential(this->liquid_potential_component);
    dealii::FEValuesExtractors::Scalar const temperature(this->temperature_component);
    dealii::QGauss<dim> quadrature_rule(fe.degree+1); // TODO: map to component...
    dealii::FEValues<dim> fe_values(fe, quadrature_rule, 
        dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);
    unsigned int const dofs_per_cell = fe.dofs_per_cell;
    unsigned int const n_q_points = quadrature_rule.size();
    dealii::Vector<double> cell_load_vector(dofs_per_cell);
    std::vector<double> solid_phase_diffusion_coefficient_values(n_q_points);
    std::vector<double> liquid_phase_diffusion_coefficient_values(n_q_points);
    std::vector<dealii::Tensor<1, dim> > solid_potential_gradients(n_q_points);
    std::vector<dealii::Tensor<1, dim> > liquid_potential_gradients(n_q_points);
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
unsigned int const n_components = dealii::DoFTools::n_components(this->dof_handler);
dealii::ComponentMask mask(n_components, false);
mask.set(this->temperature_component, true);
DoFExtractor dof_extractor(mask, mask, dofs_per_cell);
    typename dealii::DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        end_cell = this->dof_handler.end();
    for ( ; cell != end_cell; ++cell) {
        cell_load_vector = 0.0;
        fe_values.reinit(cell);
        (this->mp_values)->get_values("solid_electrical_conductivity", cell, solid_phase_diffusion_coefficient_values);
        (this->mp_values)->get_values("liquid_electrical_conductivity", cell, liquid_phase_diffusion_coefficient_values);
        fe_values[solid_potential].get_function_gradients(potential_solution_vector, solid_potential_gradients);
        fe_values[liquid_potential].get_function_gradients(potential_solution_vector, liquid_potential_gradients);
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                cell_load_vector(i) += (
                    + solid_phase_diffusion_coefficient_values[q_point] *
                    ( solid_potential_gradients[q_point] *
                      solid_potential_gradients[q_point] )
                    + liquid_phase_diffusion_coefficient_values[q_point] *
                    ( liquid_potential_gradients[q_point] *
                      liquid_potential_gradients[q_point] )
                    + coeff *
                    ( solid_phase_diffusion_coefficient_values[q_point] * 
                      solid_potential_gradients[q_point] 
                    +
                      liquid_phase_diffusion_coefficient_values[q_point] * 
                      liquid_potential_gradients[q_point] ).norm()
                    ) * fe_values[temperature].value(i, q_point) * fe_values.JxW(q_point);
            } // end for i
        } // end for quadrature point
        cell->get_dof_indices(local_dof_indices);
std::vector<dealii::types::global_dof_index> tmp_indices = dof_extractor.extract_row_indices(local_dof_indices);        
dealii::Vector<double> tmp_load_vector = dof_extractor.extract_vector(cell_load_vector);
// TODO:
//std::transform(tmp_indices.begin(), tmp_indices.end(), tmp_indices.begin(), std::bind2nd(std::minus<dealii::types::global_dof_index>(), 672*2));
        this->constraint_matrix.distribute_local_to_global(tmp_load_vector, tmp_indices, thermal_load_vector);
    } // end for cell
}

} // end namespace cache
