#define BOOST_TEST_MODULE TestDiscreteFourierTransform
#define BOOST_TEST_MAIN
#include <cap/energy_storage_device.h>
#include <cap/electrochemical_impedance_spectroscopy.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <complex>

BOOST_AUTO_TEST_CASE( test_discrete_fourier_transform )
{
    size_t const n = 32;
    size_t const k = 3;
    std::vector<double> x(n);
    for (size_t i = 0; i < n; ++i)
        x[i] = static_cast<double>(i)/(n-1);
    double const pi = acos(-1.0);
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i)
        y[i] = 1.0 + std::sin(2*pi*k*x[i]);
    auto fft_data = cap::compute_fft(&(y[0]), n);
    auto fft_freq = cap::compute_fft_frequencies(n, 1.0/n);
    for (size_t i = 0; i < fft_data.size(); ++i)
        std::cout<<i<<"  "<<fft_freq[i]<<"  "<<fft_data[i]<<"\n";
// import numpy
// n=32
// k=3
// x=numpy.linspace(0,1,n)
// y=1+numpy.sin(k*2*numpy.pi*x)
// fft_data=numpy.fft.rfft(y)
// fft_freq=numpy.fft.rfftfreq(n)
// #for i in range(len(fft_data)):
// #    print i,fft_freq[i],fft_data[i]
// print numpy.real(fft_data)
// print numpy.imag(fft_data)
// print fft_freq
    std::vector<double> real = {
        32.        ,  0.03429769,  0.21085557,  4.49628784, -0.73607911,
        -0.47868609, -0.40249877, -0.3674597 , -0.34801027, -0.3360312 ,
        -0.32817531, -0.32283154, -0.31914224, -0.3166201 , -0.31497884,
        -0.31405192, -0.31375202,
    };
    std::vector<double> imag = {
        0.        ,  -0.34823028,  -1.06004251, -14.82227458,   1.77705217,
        0.89555868,   0.60238199,   0.44775094,   0.34801027,   0.27577368,
        0.21927973,   0.17255705,   0.13219304,   0.09604566,   0.06265319,
        0.03093141,   0.        ,
    };
    std::vector<double> freq = {
        0.  ,    0.03125, 0.0625,  0.09375, 0.125,   0.15625, 0.1875,  0.21875,
        0.25,    0.28125, 0.3125,  0.34375, 0.375,   0.40625, 0.4375,  0.46875,
        0.5 ,
    };
    double const percent_tolerance = 1.0e-6;
    for (size_t i = 0; i < n/2+1; ++i)
    {
        BOOST_TEST( std::real(fft_data[i]) == real[i], boost::test_tools::tolerance(percent_tolerance) );
        BOOST_TEST( std::imag(fft_data[i]) == imag[i], boost::test_tools::tolerance(percent_tolerance) );
        BOOST_TEST( fft_freq[i]            == freq[i], boost::test_tools::tolerance(percent_tolerance) );
    }
}
