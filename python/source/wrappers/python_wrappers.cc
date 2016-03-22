/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/version.h>
#include <boost/python.hpp>

namespace pycap
{
void export_property_tree();
void export_energy_storage_device();
}

char const * pycap_docstring =
  "                                                                         \n"
  "PyCap                                                                    \n"
  "=====                                                                    \n"
  "Provides                                                                 \n"
  "  1. Energy storage device models                                        \n"
  "      Cap comes with a finite element model for supercapacitors as well  \n"
  "      as equivalents circuits models.                                    \n"
  "  2. Electrochemical measurement techniques                              \n"
  "      Cap can run a wide range of techniques on the available energy     \n"
  "      storage devices, from the simple recording of discharge curves to  \n"
  "      the more complex calculation of impedance spectra.                 \n"
  "                                                                         \n"
  "How to use the documentation                                             \n"
  "----------------------------                                             \n"
  "Documentation is available in two forms: docstrings provided with the    \n"
  "code, and a loose standing reference guide, available from               \n"
  "`Cap's online documentation <https://cap.readthedocs.org>`_.             \n"
  "                                                                         \n"
  "Use the built-in ``help`` function to view a function's docstring:       \n"
  "  >>> import pycap                                                       \n"
  "  >>> help(pycap.PropertyTree)                                           \n"
  "  ...                                                                    \n"
  "                                                                         \n"
  "Other                                                                    \n"
  "-----                                                                    \n"
  "PropertyTree                                                             \n"
  "    Wrappers for Boost.PropertyTree                                      \n"
  "    A tree data structure that stores configuration data.                \n"
  "                                                                         \n"
  "Available energy storage devices                                         \n"
  "--------------------------------                                         \n"
  "EnergyStorageDevice                                                      \n"
  "    Wrappers for Cap.EnergyStorageDevice                                 \n"
  "    See documentation.                                                   \n"
  "                                                                         \n"
  "Available electrochemical techniques                                     \n"
  "------------------------------------                                     \n"
  "Discharge                                                                \n"
  "    Performs the discharge in one of four different control modes.       \n"
  "Charge                                                                   \n"
  "    Performs the charge with an optional constant voltage step at the    \n"
  "    end.                                                                 \n"
  "CyclicChargeDischarge                                                    \n"
  "    Records charge and discharge curves through a number of cycles.      \n"
  "CyclicVoltammetry                                                        \n"
  "    Applies cyclic linear voltage ramps.                                 \n"
  "ElectrochemicalImpedanceSpectroscopy                                     \n"
  "    Measures the complex electrochemaical impedance over a range of      \n"
  "    frequencies.                                                         \n"
  "RagoneAnalysis                                                           \n"
  "    Complete a series of discharges at various rate to produce a ragone  \n"
  "    plot.                                                                \n"
  ;

BOOST_PYTHON_MODULE(PyCap)
{
  boost::python::scope().attr("__version__"        ) = cap::version()        ;
  boost::python::scope().attr("__git_branch__"     ) = cap::git_branch()     ;
  boost::python::scope().attr("__git_commit_hash__") = cap::git_commit_hash();
  boost::python::scope().attr("__git_remote_url__" ) = cap::git_remote_url() ;
  boost::python::scope().attr("__doc__") = pycap_docstring;

  boost::python::docstring_options doc_options;
  doc_options.enable_user_defined();
  doc_options.enable_py_signatures();
  doc_options.disable_cpp_signatures();

  pycap::export_energy_storage_device();

  pycap::export_property_tree();
}

