#define BOOST_TEST_MODULE ExactTransientSolution
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/mp_values.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

namespace cap {

void compute_parameters(std::shared_ptr<boost::property_tree::ptree const> input_database,
                        std::shared_ptr<boost::property_tree::ptree      > output_database)
{
    double const sandwich_height = input_database->get<double>("geometry.sandwich_height");
    double const cross_sectional_area = sandwich_height * 1.0;
    double const electrode_width = input_database->get<double>("geometry.electrode_width");
    double const separator_width = input_database->get<double>("geometry.separator_width");

    // getting the material parameters values
    std::shared_ptr<boost::property_tree::ptree> material_properties_database = 
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("material_properties"));
    std::shared_ptr<cap::SuperCapacitorMPValues<2> > mp_values = std::shared_ptr<cap::SuperCapacitorMPValues<2> >
        (new cap::SuperCapacitorMPValues<2>(cap::SuperCapacitorMPValuesParameters<2>(material_properties_database)));
    // build dummy cell itertor and set its material id
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube (triangulation);
    dealii::DoFHandler<2> dof_handler(triangulation);
    dealii::DoFHandler<2>::active_cell_iterator cell = 
        dof_handler.begin_active();
    // electrode
    cell->set_material_id(
        input_database->get<dealii::types::material_id>("material_properties.anode_electrode_material_id"));
    std::vector<double> electrode_solid_electrical_conductivity_values(1);
    std::vector<double> electrode_liquid_electrical_conductivity_values(1);
    std::vector<double> electrode_specific_capacitance_values(1);
    std::vector<double> electrode_exchange_current_density_values(1);
    std::vector<double> electrode_electron_thermal_voltage_values(1);
    mp_values->get_values("solid_electrical_conductivity" , cell, electrode_solid_electrical_conductivity_values );
    mp_values->get_values("liquid_electrical_conductivity", cell, electrode_liquid_electrical_conductivity_values);
    mp_values->get_values("specific_capacitance"          , cell, electrode_specific_capacitance_values          );
    mp_values->get_values("faradaic_reaction_coefficient" , cell, electrode_exchange_current_density_values      );
    mp_values->get_values("electron_thermal_voltage"      , cell, electrode_electron_thermal_voltage_values      );
    double const total_current = -1.0; // normalized
    double const dimensionless_exchange_current_density = electrode_exchange_current_density_values[0]
        * std::pow(electrode_width,2) 
        * ( 1.0 / electrode_solid_electrical_conductivity_values[0]
          + 1.0 / electrode_liquid_electrical_conductivity_values[0] );
    double const dimensionless_current_density = total_current
        * electrode_width
        / electrode_liquid_electrical_conductivity_values[0]
        / electrode_electron_thermal_voltage_values[0];
    double const ratio_of_solution_phase_to_matrix_phase_conductivities =
        electrode_liquid_electrical_conductivity_values[0] / electrode_solid_electrical_conductivity_values[0];

    output_database->put("dimensionless_current_density"                         , dimensionless_current_density                         );
    output_database->put("dimensionless_exchange_current_density"                , dimensionless_exchange_current_density                );
    output_database->put("ratio_of_solution_phase_to_matrix_phase_conductivities", ratio_of_solution_phase_to_matrix_phase_conductivities);

    output_database->put("position_normalization_factor", electrode_width);
    output_database->put("time_normalization_factor"    ,
        electrode_specific_capacitance_values[0] 
            * ( 1.0 / electrode_solid_electrical_conductivity_values[0]
              + 1.0 / electrode_liquid_electrical_conductivity_values[0] )
            * std::pow(electrode_width,2) );
          
    // separator
    cell->set_material_id(
        input_database->get<dealii::types::material_id>("material_properties.separator_material_id"));
    std::vector<double> separator_liquid_electrical_conductivity_values(1);
    mp_values->get_values("liquid_electrical_conductivity", cell, separator_liquid_electrical_conductivity_values);

    double const potential_drop_across_the_separator = -total_current * separator_width / separator_liquid_electrical_conductivity_values[0];
    double const voltage_normalization_factor        = electrode_electron_thermal_voltage_values[0];
    output_database->put("potential_drop_across_the_separator", potential_drop_across_the_separator);
    output_database->put("voltage_normalization_factor"       , voltage_normalization_factor       );
    output_database->put("cross_sectional_area"               , cross_sectional_area               );
}

void verification_problem(std::shared_ptr<cap::EnergyStorageDevice> dev, std::shared_ptr<boost::property_tree::ptree const> database, std::ostream & os = std::cout)
{
    double dimensionless_current_density                          = database->get<double>("dimensionless_current_density"                         );
    double dimensionless_exchange_current_density                 = database->get<double>("dimensionless_exchange_current_density"                );
    double ratio_of_solution_phase_to_matrix_phase_conductivities = database->get<double>("ratio_of_solution_phase_to_matrix_phase_conductivities");

    int    const infty = database->get<int>("terms_in_truncation_of_infinite_series");
    double const pi    = std::acos(-1.0);

//    auto compute_dimensionless_time =
//        [&electrode_solid_electrical_conductivity_values, &electrode_liquid_electrical_conductivity_values,
//        &electrode_specific_capacitance_values, electrode_width]
//        (double const time)
//        {
//            return time / 
//            ( electrode_specific_capacitance_values[0] 
//            * ( 1.0 / electrode_solid_electrical_conductivity_values[0]
//              + 1.0 / electrode_liquid_electrical_conductivity_values[0] )
//            * std::pow(electrode_width,2) );
//        };
//    auto compute_dimensionless_position =
//        [electrode_width](double const position) { return position / electrode_width; };

    auto compute_dimensionless_overpotential =
        [infty, pi,
        &ratio_of_solution_phase_to_matrix_phase_conductivities, &dimensionless_exchange_current_density, &dimensionless_current_density]
        (double const dimensionless_time, double const dimensionless_position)
        {
            std::vector<double> coefficients(infty);
            for (int n = 0; n < infty; ++n) {
                coefficients[n] =
                ( ratio_of_solution_phase_to_matrix_phase_conductivities * std::cos( n * pi ) + 1.0 )
                / ( dimensionless_exchange_current_density + std::pow(n,2) * std::pow(pi,2) )
                * std::cos( n * pi * dimensionless_position )
                * std::exp( - (std::pow(n,2) * std::pow(pi,2) + dimensionless_exchange_current_density) * dimensionless_time );
            }
            return 
                dimensionless_current_density 
                    * ( 1.0 + ratio_of_solution_phase_to_matrix_phase_conductivities )
                    * std::exp( - dimensionless_exchange_current_density * dimensionless_time )
                    / dimensionless_exchange_current_density
                -
                dimensionless_current_density
                    * ( std::cosh( std::sqrt(dimensionless_exchange_current_density) * (1.0 - dimensionless_position) )
                      +
                        ratio_of_solution_phase_to_matrix_phase_conductivities * std::cosh( std::sqrt(dimensionless_exchange_current_density) * dimensionless_position ) )
                    / ( std::sqrt(dimensionless_exchange_current_density) * std::sinh( std::sqrt(dimensionless_exchange_current_density) ) )
                +
                2.0 * dimensionless_current_density * std::accumulate(&(coefficients[1]), &(coefficients[infty]), 0.0);
        };

    auto compute_dimensionless_potential_drop_across_the_electrode =
        [&compute_dimensionless_overpotential,
        &ratio_of_solution_phase_to_matrix_phase_conductivities, &dimensionless_current_density]
        (double const dimensionless_time)
        {
            return
            ( compute_dimensionless_overpotential(dimensionless_time, 0.0)
            + ratio_of_solution_phase_to_matrix_phase_conductivities * compute_dimensionless_overpotential(dimensionless_time, 1.0)
            - dimensionless_current_density * ratio_of_solution_phase_to_matrix_phase_conductivities
            ) / ( 1.0 + ratio_of_solution_phase_to_matrix_phase_conductivities);
        };

    auto compute_dimensionless_interfacial_current_density =
        [infty, pi,
        &ratio_of_solution_phase_to_matrix_phase_conductivities, &dimensionless_exchange_current_density, &dimensionless_current_density]
        (double const dimensionless_time, double const dimensionless_position)
        {
            std::vector<double> coefficients(infty);
            for (int n = 0; n < infty; ++n) {
                coefficients[n] =
                ( ratio_of_solution_phase_to_matrix_phase_conductivities * std::cos( n * pi ) + 1.0 )
                / ( dimensionless_exchange_current_density + std::pow(n,2) * std::pow(pi,2) )
                / (ratio_of_solution_phase_to_matrix_phase_conductivities + 1.0)
                * std::pow(n,2) * std::pow(pi,2) * std::cos( n * pi * dimensionless_position )
                * std::exp( - (std::pow(n,2) * std::pow(pi,2) + dimensionless_exchange_current_density) * dimensionless_time );
            }
            return
                - 2.0 * dimensionless_current_density * std::accumulate(&(coefficients[1]), &(coefficients[infty]), 0.0)
                - dimensionless_current_density * std::sqrt(dimensionless_exchange_current_density)
                    * ( std::cosh( std::sqrt(dimensionless_exchange_current_density) * (1.0 - dimensionless_position) )
                      +
                        ratio_of_solution_phase_to_matrix_phase_conductivities * std::cosh( std::sqrt(dimensionless_exchange_current_density) * dimensionless_position ) )
                    / ( (ratio_of_solution_phase_to_matrix_phase_conductivities + 1.0) * std::sinh( std::sqrt(dimensionless_exchange_current_density) ) );
        };


    // exact vs computed
    double const charge_current = database->get<double>("charge_current");
    double const charge_time    = database->get<double>("charge_time"   );
    double const time_step      = database->get<double>("time_step"     );
    double const epsilon        = time_step * 1.0e-4;
    double const cross_sectional_area                = database->get<double>("cross_sectional_area"               );
    double const time_normalization_factor           = database->get<double>("time_normalization_factor"          );
    double const voltage_normalization_factor        = database->get<double>("voltage_normalization_factor"       );
    double       potential_drop_across_the_separator = database->get<double>("potential_drop_across_the_separator");
    potential_drop_across_the_separator *= charge_current / cross_sectional_area;
    dimensionless_current_density *= charge_current / cross_sectional_area;

    std::cout<<"delta="<<dimensionless_current_density<<"\n";
    std::cout<<"nu2="<<dimensionless_exchange_current_density<<"\n";
    std::cout<<"beta="<<ratio_of_solution_phase_to_matrix_phase_conductivities<<"\n";

    dev->reset_voltage(0.0);
    double computed_voltage;
    double exact_voltage;
    for (double time = 0.0; time <= charge_time+epsilon; time += time_step)
    {
        double const dimensionless_time = (time+time_step) / time_normalization_factor;
        double const dimensionless_potential_drop_across_the_electrode = compute_dimensionless_potential_drop_across_the_electrode(dimensionless_time);
        exact_voltage = 2.0 * dimensionless_potential_drop_across_the_electrode * voltage_normalization_factor + potential_drop_across_the_separator;
        dev->evolve_one_time_step_constant_current(time_step, charge_current);
        dev->get_voltage(computed_voltage);
        os<<boost::format("  %22.15e  %22.15e  %22.15e  \n")
                % time 
                % exact_voltage
                % computed_voltage
                ;
    }
    double const percent_tolerance = database->get<double>("percent_tolerance");
    BOOST_CHECK_CLOSE(computed_voltage, exact_voltage, percent_tolerance);

    // figure 2
    dimensionless_current_density                          = -1.0; 
    dimensionless_exchange_current_density                 = 1.0;
    ratio_of_solution_phase_to_matrix_phase_conductivities = 1.0;
    std::vector<double> tau { 0.1, 0.5, 1.0, 2.5 };
    int const n = 100;
    std::vector<double> x(n+1);
    for (int i = 0; i < n+1; ++i)
        x[i] = static_cast<double>(i) / n;
    std::fstream fout;
    for (double const & dimensionless_time : tau)
    {
        fout.open("fig2_"+std::to_string(dimensionless_time), std::fstream::out);
        for (double const & dimensionless_position : x) 
        {
            double const dimensionless_overpotential = compute_dimensionless_overpotential(dimensionless_time, dimensionless_position);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_position
                % dimensionless_overpotential
                ;
        }
        fout.close();
    }

    // figure 3
    tau.clear();
    tau.resize(n+1);
    for (int i = 0; i < n+1; ++i)
        tau[i] = 10.0 * static_cast<double>(i) / n;

    for (double const & delta : { 0.1, 0.5, 1.0, 2.5 })
    {
        dimensionless_current_density                          = -delta;
        dimensionless_exchange_current_density                 = 1.0;
        ratio_of_solution_phase_to_matrix_phase_conductivities = 1.0;
        fout.open("fig3_"+std::to_string(delta), std::fstream::out);
        for (double const & dimensionless_time : tau)
        {
            double const dimensionless_potential_drop = compute_dimensionless_potential_drop_across_the_electrode(dimensionless_time);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_time
                % dimensionless_potential_drop
                ;
        }
        fout.close();
    }

    // figure 4
    tau.clear();
    tau.resize(n+1);
    for (int i = 0; i < n+1; ++i)
        tau[i] = 5.0 * static_cast<double>(i) / n;
    for (double const & delta : { 1.0, 2.5 })
    {
        for (double const & nu2 : { 0.01, 1.0 })
        {
            dimensionless_current_density                          = -delta;
            dimensionless_exchange_current_density                 = nu2;
            ratio_of_solution_phase_to_matrix_phase_conductivities = 0.0;
            fout.open("fig4_"+std::to_string(delta)+"_"+std::to_string(nu2), std::fstream::out);
            for (double const & dimensionless_time : tau)
            {
                double const dimensionless_potential_drop = compute_dimensionless_potential_drop_across_the_electrode(dimensionless_time);
                fout<<boost::format("  %22.15e  %22.15e  \n")
                    % dimensionless_time
                    % dimensionless_potential_drop
                    ;
            }
            fout.close();
        }
    }

    // figure 5
    dimensionless_current_density                          = -1.0;
    dimensionless_exchange_current_density                 = 1.0;
    ratio_of_solution_phase_to_matrix_phase_conductivities = 1.0;
    tau.clear();
    tau = std::vector<double> { 0.01, 0.02, 0.04, 0.1, 0.15, 0.2 };
    for (double const & dimensionless_time : tau)
    {
        fout.open("fig5_"+std::to_string(dimensionless_time), std::fstream::out);
        for (double const & dimensionless_position : x) 
        {
            double const dimensionless_interfacial_current_density = compute_dimensionless_interfacial_current_density(dimensionless_time, dimensionless_position);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_position
                % dimensionless_interfacial_current_density
                ;
        }
        fout.close();
    }

    // figure 6
    dimensionless_current_density                          = -1.0;
    dimensionless_exchange_current_density                 = 0.01;
    ratio_of_solution_phase_to_matrix_phase_conductivities = 1.0;
    for (double const & dimensionless_time : tau)
    {
        fout.open("fig6_"+std::to_string(dimensionless_time), std::fstream::out);
        for (double const & dimensionless_position : x) 
        {
            double const dimensionless_interfacial_current_density = compute_dimensionless_interfacial_current_density(dimensionless_time, dimensionless_position);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_position
                % dimensionless_interfacial_current_density
                ;
        }
        fout.close();
    }

    // figure 7
    for (double const & delta : { 0.5, 1.0, 2.0, 2.5 })
    {
        dimensionless_current_density                          = -delta;
        dimensionless_exchange_current_density                 = 1.0;
        ratio_of_solution_phase_to_matrix_phase_conductivities = 1.0;
        fout.open("fig7_"+std::to_string(delta), std::fstream::out);
        double const dimensionless_time = 0.2;
        for (double const & dimensionless_position : x) 
        {
            double const dimensionless_interfacial_current_density = compute_dimensionless_interfacial_current_density(dimensionless_time, dimensionless_position);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_position
                % dimensionless_interfacial_current_density
                ;
        }
        fout.close();
    }

    // figure 8
    dimensionless_current_density                          = -1.0;
    dimensionless_exchange_current_density                 = 1.0;
    ratio_of_solution_phase_to_matrix_phase_conductivities = 0.0;
    for (double const & dimensionless_time : tau)
    {
        fout.open("fig8_"+std::to_string(dimensionless_time), std::fstream::out);
        for (double const & dimensionless_position : x) 
        {
            double const dimensionless_interfacial_current_density = compute_dimensionless_interfacial_current_density(dimensionless_time, dimensionless_position);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_position
                % dimensionless_interfacial_current_density
                ;
        }
        fout.close();
    }

    // figure 9
    dimensionless_current_density                          = -1.0;
    dimensionless_exchange_current_density                 = 1.0;
    ratio_of_solution_phase_to_matrix_phase_conductivities = 10.0;
    for (double const & dimensionless_time : tau)
    {
        fout.open("fig9_"+std::to_string(dimensionless_time), std::fstream::out);
        for (double const & dimensionless_position : x) 
        {
            double const dimensionless_interfacial_current_density = compute_dimensionless_interfacial_current_density(dimensionless_time, dimensionless_position);
            fout<<boost::format("  %22.15e  %22.15e  \n")
                % dimensionless_position
                % dimensionless_interfacial_current_density
                ;
        }
        fout.close();
    }
}

} // end namespace cap

BOOST_AUTO_TEST_CASE( test_exact_transient_solution )
{
    // parse input file
    std::shared_ptr<boost::property_tree::ptree> input_database =
        std::make_shared<boost::property_tree::ptree>();
    read_xml("input_verification_problem", *input_database);

    // build an energy storage system
    std::shared_ptr<boost::property_tree::ptree> device_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("device"));
    std::shared_ptr<cap::EnergyStorageDevice> device =
        cap::buildEnergyStorageDevice(std::make_shared<cap::Parameters>(device_database));

    // measure discharge curve
    std::fstream fout;
    fout.open("verification_problem_data", std::fstream::out);

    std::shared_ptr<boost::property_tree::ptree> verification_problem_database =
        std::make_shared<boost::property_tree::ptree>(input_database->get_child("verification_problem"));

    cap::compute_parameters(device_database, verification_problem_database);

    cap::verification_problem(device, verification_problem_database, fout);

    fout.close();
}    
