/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_THERMAL_OPERATOR_H
#define CAP_THERMAL_OPERATOR_H

#include <cap/operator.h>

namespace cap
{

//////////////////////// THERMAL OPERATOR PARAMETERS ////////////////////////
/**
 * This class encapsulates the parameters used ThermalOperator.
 */
template <int dim>
class ThermalOperatorParameters : public OperatorParameters<dim>
{
public:
  [[deprecated]] ThermalOperatorParameters(
      std::shared_ptr<boost::property_tree::ptree const> d)
      : OperatorParameters<dim>(d)
  {
  }
};

//////////////////////// THERMAL OPERATOR ////////////////////////////
/**
 * Create the system of equations associated to the heat transfer. The heat
 * source is computed in ElectrochemicalOperator.
 */
template <int dim>
class ThermalOperator : public Operator<dim>
{
public:
  [[deprecated]] ThermalOperator(
      std::shared_ptr<OperatorParameters<dim> const> parameters);

  void
  reset(std::shared_ptr<OperatorParameters<dim> const> parameters) override;

protected:
  void compute_thermal_operator_contribution();
  void compute_robin_boundary_contribution();

  unsigned int temperature_component;
};

} // end namespace cap

#endif // CAP_THERMAL_OPERATOR_H
