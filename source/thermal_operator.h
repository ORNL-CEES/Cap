#ifndef CACHE_THERMAL_OPERATOR_H
#define CACHE_THERMAL_OPERATOR_H

#include <cache/operator.h>

namespace cache {

//////////////////////// THERMAL OPERATOR PARAMETERS ////////////////////////////
template <int dim>
class ThermalOperatorParameters : public OperatorParameters<dim> {
public:
    unsigned int temperature_component;
};

//////////////////////// THERMAL OPERATOR ////////////////////////////

template <int dim>
class ThermalOperator : public Operator<dim> {
public:
    ThermalOperator(OperatorParameters<dim> const & parameters);
    void reset(OperatorParameters<dim> const & parameters);
    
protected:
    void compute_thermal_operator_contribution();
    void compute_robin_boundary_contribution();

    unsigned int temperature_component;
};

} // end namespace cache

#endif // CACHE_THERMAL_OPERATOR_H
