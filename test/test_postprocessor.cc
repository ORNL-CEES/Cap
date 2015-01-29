#include <boost/format.hpp>
#include <vector>
#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cmath>

namespace cap {

template <class InputIterator1, class InputIterator2,
          class OutputIterator, class T>
void approximate_integral_with_trapezoidal_rule
    (InputIterator1 first1, InputIterator1 last1,
                            InputIterator2 first2,
     OutputIterator result, T init)
{
    *result = init;
    ++first1; ++first2; ++result;
    while (first1 != last1)
    {
        *result = *std::prev(result) + 0.5 * (*first1 - *std::prev(first1)) * (*first2 + *std::prev(first2));
        ++first1; ++first2; ++result;
    }

}

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

} // end namespace cap

//  #!/usr/bin/env python
//  from pylab import *
//  data=loadtxt("data",usecols=(0,2,3))
//  time=data[:,0]
//  power=data[:,1]
//  energy=data[:,2]
//  fig,axarr=subplots(2,sharex=True)
//  axarr[0].plot(time,power)
//  axarr[0].set_ylabel("power")
//  axarr[1].plot(time,energy)
//  axarr[1].set_xlabel("time")
//  axarr[1].set_ylabel("energy")
//  show()

int main(int argc, char *argv[])
{
    std::fstream fout("data", std::fstream::out);
    std::ostream & os = fout;
    bool const verbose = true;

    double const initial_time = 0.0;
    double const final_time   = 2.0;
    std::size_t const n = 10001;
    double const pi = 3.14159265359;
    std::vector<double> time(n);
    std::vector<double> power(n);
    std::vector<std::string> capacitor_state(n);

    // INITIALIZE SOLUTION
    for (std::size_t i = 0; i < n; ++i) {
        time[i] = initial_time + static_cast<double>(i) / (n - 1) * (final_time - initial_time);
        power[i] = std::cos(2.0 * pi * time[i]);
        capacitor_state[i] = (power[i] > 0.0 ? "charging" : "discharging");
    }

    // COMPUTE ENERGY
    std::vector<double> energy(n); 
    cap::compute_energy(capacitor_state, time, power, energy);

    // GET DURATION AND AVERAGED POWER
    std::vector<double> duration;
    std::vector<double> average_power;
    cap::extract_duration_and_average_power(capacitor_state, time, energy, duration, average_power);

    for (std::size_t i = 0; i < duration.size(); ++i)
        std::cout<<duration[i]<<"  "<<average_power[i]<<"\n";

    // CHECK THE ANSWER
    std::vector<double> exact(n);
    std::vector<double> error(n);
    double correction = 0.0;
    exact[0] = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        exact[i] = 0.5 / pi * std::sin(2.0 * pi * time[i]) - correction;
        error[i] = energy[i] - exact[i];
    }

    // PRINT TO THE STREAM
    if (verbose)
        for (std::size_t i = 0; i < n; ++i)
            os<<boost::format("  %10.5f  %15s  %10.5f  %10.5f  %10.5f  %10.5f \n")
                % time[i]
                % capacitor_state[i]
                % power[i]
                % energy[i]
                % exact[i]
                % error[i]
                ;
    
    return 0;
}
