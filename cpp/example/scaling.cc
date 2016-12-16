#include <cap/energy_storage_device.h>
#include <cap/mp_values.h>
#include <cap/default_inspector.h>
#include <cap/supercapacitor.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/timer.hpp>
#include <iostream>

void run_example(boost::mpi::communicator &comm)
{
  boost::mpi::timer timer;
  if (comm.rank() == 0)
    std::cout << "Number of processors: " << comm.size() << std::endl;

  // Parse input file
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database, comm);

  unsigned int const n_time_steps = 10;
  double const time_step = 0.1;
  double const charge_voltage = 2.1;
  for (unsigned int i = 0; i < n_time_steps; ++i)
    device->evolve_one_time_step_constant_voltage(time_step, charge_voltage);

  if (comm.rank() == 0)
  {
    cap::DefaultInspector inspector;
    inspector.inspect(device.get());
    auto data = inspector.get_data();
    std::cout << "n dofs: " << data["n_dofs"] << std::endl;
    std::cout << "Elapsed time: " << timer.elapsed() << std::endl;
  }

  if (device_database.get<int>("dim") == 2)
  {
    cap::SuperCapacitorInspector<2> supercap_inspector;
    supercap_inspector.inspect(device.get());
  }
  else
  {
    cap::SuperCapacitorInspector<3> supercap_inspector;
    supercap_inspector.inspect(device.get());
  }
}

int main(int argc, char *argv[])
{
  try
  {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    run_example(world);
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
