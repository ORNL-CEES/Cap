#ifndef CACHE_THERMAL_OPERATOR_H
#define CACHE_THERMAL_OPERATOR_H

#include <cache/operator.h>

namespace cache {

//////////////////////// THERMAL OPERATOR PARAMETERS ////////////////////////////
template <int dim>
class ThermalOperatorParameters : public OperatorParameters<dim> {
public:
    ThermalOperatorParameters(std::shared_ptr<boost::property_tree::ptree const> d)
        : OperatorParameters<dim>(d)
    { }
};

//////////////////////// THERMAL OPERATOR ////////////////////////////

template <int dim>
class ThermalOperator : public Operator<dim> {
public:
    ThermalOperator(std::shared_ptr<OperatorParameters<dim> const> parameters);
    void reset(std::shared_ptr<OperatorParameters<dim> const> parameters);
    
protected:
    void compute_thermal_operator_contribution();
    void compute_robin_boundary_contribution();

    unsigned int temperature_component;
};

} // end namespace cache

#endif // CACHE_THERMAL_OPERATOR_H
