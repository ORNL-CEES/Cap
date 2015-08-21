#include <cap/electrochemical_impedance_spectroscopy.h>
#include <cap/utils.h>
#include <boost/math/special_functions/sin_pi.hpp>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#include <iostream>

namespace cap {

std::vector<std::complex<double>>
compute_fft(double const * input, size_t const n)
{
    std::vector<double> data(input, input+n);
    gsl_fft_real_radix2_transform(&(data[0]), 1, n);
    std::vector<double> unpacked_data(2*n);
    gsl_fft_halfcomplex_radix2_unpack(&(data[0]), &(unpacked_data[0]), 1, n);
    std::vector<std::complex<double>> output(n/2+1);
    for (size_t k = 0; k < n/2+1; ++k)            
        output[k] = std::complex<double>(REAL(unpacked_data,k), IMAG(unpacked_data,k));
    return output;
}

std::vector<double>
compute_fft_frequencies(size_t const n, double const frequency)
{
    std::vector<double> fft_frequencies(n/2+1);
    for (size_t k = 0; k < n/2+1; ++k)
        fft_frequencies[k] = k * frequency;
    return fft_frequencies;
}

std::function<void(std::shared_ptr<cap::EnergyStorageDevice>,double,double)>
get_evolve_one_time_step(std::shared_ptr<boost::property_tree::ptree const> database)
{
    std::vector<int   > const harmonics  = cap::to_vector<int   >(database->get<std::string>("harmonics" ));
    std::vector<double> const amplitudes = cap::to_vector<double>(database->get<std::string>("amplitudes"));
    std::vector<double> const phases     = cap::to_vector<double>(database->get<std::string>("phases"    ));
    double const frequency = database->get<double>("frequency");
    auto compute_ac_excitation_signal =
        [frequency,harmonics,amplitudes,phases](double time)
        {
            double excitation_signal = 0.0;
            for (size_t k = 0; k < harmonics.size(); ++k)
                excitation_signal += amplitudes[k] * boost::math::sin_pi(2 * harmonics[k] * frequency * time + phases[k]);
            return excitation_signal;
        };

    std::string const mode = database->get<std::string>("mode");
    if (mode.compare("galvanostatic") != 0) {
        return 
            [compute_ac_excitation_signal]
            (std::shared_ptr<cap::EnergyStorageDevice> device, double time, double time_step)
            {
                device->evolve_one_time_step_changing_current(time_step, compute_ac_excitation_signal(time));
            };
    } else if (mode.compare("potentiostatic") != 0) {
        double const dc_voltage = database->get<double>("dc_voltage");
        return 
            [dc_voltage,compute_ac_excitation_signal]
            (std::shared_ptr<cap::EnergyStorageDevice> device, double time, double time_step)
            {
                device->evolve_one_time_step_changing_voltage(time_step, dc_voltage+compute_ac_excitation_signal(time));
            };
    } else {
        throw std::runtime_error("invalide mode "+mode+" for EIS measurement");
    }
}

std::map<double,std::complex<double>>
measure_impedance(std::shared_ptr<cap::EnergyStorageDevice> device, std::shared_ptr<boost::property_tree::ptree const> database)
{
    std::vector<int   > const harmonics  = cap::to_vector<int   >(database->get<std::string>("harmonics" ));
    double const frequency       = database->get<double>("frequency"    );
    int    const cycles          = database->get<int   >("cycles"         );
    int    const ignore_cycles   = database->get<int   >("ignore_cycles"  );
    int    const steps_per_cycle = database->get<int   >("steps_per_cycle");
    auto evolve_one_time_step = get_evolve_one_time_step(database);

    // apply excitation signal and measure response
    std::vector<double> time(cycles*steps_per_cycle);
    std::vector<double> current(cycles*steps_per_cycle);
    std::vector<double> voltage(cycles*steps_per_cycle);
    double const time_step = 1.0 / frequency / steps_per_cycle;
    for (int step = 0; step < cycles*steps_per_cycle; ++step)
    {
          time[step] = (step+1) * time_step;
          evolve_one_time_step(device, time[step], time_step);
          device->get_current(current[step]);
          device->get_voltage(voltage[step]);
    }

    // compute discrete fourrier transforms
    std::vector<std::complex<double>> fft_current =
        compute_fft(&(current[ignore_cycles*steps_per_cycle]), (cycles-ignore_cycles)*steps_per_cycle);
    std::vector<std::complex<double>> fft_voltage =
        compute_fft(&(voltage[ignore_cycles*steps_per_cycle]), (cycles-ignore_cycles)*steps_per_cycle);
    std::vector<double> fft_frequencies =
        compute_fft_frequencies((cycles-ignore_cycles)*steps_per_cycle, frequency/(cycles-ignore_cycles));

    // check signal-to-noise ratio
    int const n = (cycles-ignore_cycles)*steps_per_cycle;
    assert(fft_current.size() == n/2+1);
    assert(fft_voltage.size() == n/2+1);
    assert(fft_frequencies.size() == n/2+1);
    std::vector<int> excited_harmonics;
    std::vector<int> unexcited_harmonics;
// TODO
    for (int i = 1; i < n/2; ++i)
        if ((std::find(harmonics.begin(), harmonics.end(), i/(cycles-ignore_cycles)) != harmonics.end())
            && (i%(cycles-ignore_cycles) == 0))
            excited_harmonics.push_back(i);
        else
            unexcited_harmonics.push_back(i);
     
    // return complex impedance
    std::map<double,std::complex<double>> impedance;
    for (int const & k : excited_harmonics)
    {
        impedance.emplace(
            fft_frequencies[k],
            fft_voltage[k] / fft_current[k]
            );
    }
    return impedance;
}

} // end namespace cap
