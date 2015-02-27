#ifndef CAP_EQUIVALENT_CIRCUIT_H
#define CAP_EQUIVALENT_CIRCUIT_H

#include <cap/resistor_capacitor.h>
#include <cap/mp_values.h>
#include <cap/electrochemical_operator.h> // TODO: has definition of the enum CapacitorState
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>

namespace cap {
enum OutputData { TEMPERATURE, VOLTAGE, CURRENT, JOULE_HEATING, SURFACE_AREA, VOLUME, MASS, N_DATA};
class EquivalentCircuit
{
public:
    EquivalentCircuit(std::shared_ptr<boost::property_tree::ptree const> datatabase);
    void process_solution(double * data);
    void evolve_one_time_step(double const & delta_t);
    void reset(std::shared_ptr<boost::property_tree::ptree const> database);
    void setup(CapacitorState const change);
    void run(std::shared_ptr<boost::property_tree::ptree const> input_params,
        std::shared_ptr<boost::property_tree::ptree> output_params);

private:
    std::shared_ptr<cap::SeriesRC> equivalent_circuit;
    cap::CapacitorState capacitor_state;
    double I_charge;
    double I_discharge;
    double U_charge;
    double U_discharge;
};

} // end namespace cap

#endif // CAP_EQUIVALENT_CIRCUIT_H
