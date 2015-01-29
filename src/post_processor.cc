#include <cap/post_processor.templates.h>

namespace cap {

template class PostprocessorParameters<2>;
template class Postprocessor<2>;
template class PostprocessorParameters<3>;
template class Postprocessor<3>;

template class SuperCapacitorPostprocessorParameters<2>;
template class SuperCapacitorPostprocessor<2>;
template class SuperCapacitorPostprocessorParameters<3>;
template class SuperCapacitorPostprocessor<3>;

void compute_energy(std::vector<std::string> const & capacitor_state,
    std::vector<double> const & time,
    std::vector<double> const & power,
    std::vector<double> & energy)
   
{
    bool const valid_input = 
           (time.size() == power.size())
        && (time.size() == energy.size())
        && (time.size() == capacitor_state.size())
        ;
    if (!valid_input) throw std::runtime_error("invalid input");
    std::vector<std::string>::const_iterator it     = capacitor_state.begin();
    std::vector<std::string>::const_iterator end_it = capacitor_state.end();
    std::size_t first =  0;
    while (it != end_it)
    {
        auto same = [&it] (std::string const & o) { return it->compare(o) == 0; };
        std::vector<std::string>::const_iterator next = std::find_if_not(it, end_it, same);
        std::size_t last = first + std::distance(it, next);
        approximate_integral_with_trapezoidal_rule(std::next(time.begin(), first), std::next(time.begin(), last),
                                                   std::next(power.begin(), first), 
                                                   std::next(energy.begin(), first), 0.0);
        it = next;
        first = last;
    }
}

void extract_duration_and_average_power(std::vector<std::string> const & capacitor_state,
    std::vector<double> const & time,
    std::vector<double> const & energy,
    std::vector<double> & duration,
    std::vector<double> & average_power)
{
    bool const valid_input = (time.size() == energy.size())
        && duration.empty()
        && average_power.empty()
        ;
    if (!valid_input) throw std::runtime_error("invalid input");
    std::vector<std::string>::const_iterator it     = capacitor_state.begin();
    std::vector<std::string>::const_iterator end_it = capacitor_state.end();
    std::size_t first =  0;
    while (it != end_it)
    {
        auto same = [&it] (std::string const & o) { return it->compare(o) == 0; };
        std::vector<std::string>::const_iterator next = std::find_if_not(it, end_it, same);
        std::size_t last = first + std::distance(it, next);
        duration.push_back(time[last-1] - time[first]);
        average_power.push_back((energy[last-1] - energy[first]) / duration.back());
        it = next;
        first = last;
    }
}

void compute_thermal_energy_losses(std::vector<std::string> const & capacitor_state,
    std::vector<double> const & time,
    std::vector<double> const & heat_production,
    std::vector<double> & energy_losses)
{
    bool const valid_input = (time.size() == energy_losses.size())
        && (time.size() == heat_production.size())
        && (time.size() == capacitor_state.size())
        ;
    if (!valid_input) throw std::runtime_error("invalid input");
    approximate_integral_with_trapezoidal_rule(time.begin(), time.end(),
        heat_production.begin(), energy_losses.begin(), 0.0);
}

} // end namespace cap
