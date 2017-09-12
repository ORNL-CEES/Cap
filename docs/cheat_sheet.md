Parameters 
----------
You can find below a list of the parameters that describe a supercapacitor. Not
all of them can be used at the same time and some have default values:
  1. solid_potential_component (unsigned int)
  2. liquid_potential_component (unsigned int)
  3. dim (unsigned int)
  4. geometry
    * divisions (string)
    * dimensions (string)
    * type (string)
    * mesh_file(string)
    * coarse_mesh_filename (string)
    * materials (int)
    * material_id (string)
    * coarse_mesh_filename (string)
    * materials (int)
    * material_id (string)
    * materials (int)
    * material_X (X in [0, materials))
      a. material_id (string)
      b. name (string)
      c. weight (unsigned int)
    * boundaries (int)
    * boundary_X (X in [0, boundaries))
      a. boundary_id (string)
      b. name (string)
    * anode_collector_thickness (double)
    * anode_electrode_thickness (double)
    * separator_thickness (double)
    * cathode_collector_thickness (double)
    * cathode_electrode_thickness (double)
    * geometric_area (double)
    * tab_height (double)
    * n_refinements (unsigned int)
    * shape (string)
    * checkpoint (bool)
    * n_repetitions (unsigned int)
  5. material_properties
    * material_name
      a. type (string)
      b. matrix_phase (string)
      c. solution_phase (string)
      d. metal_foil (string)
      e. anodic_charge_transfer_coefficient (double)
      f. cathodic_charge_transfer_coefficient (double)
      g. faraday_constant (double)
      h. gas_constant (double)
      i. temperature (double)
    * matrix_phase_X
      a. differential_capacitance
      b. exchange_current_density
      c. electrical_resistivity
      d. void_volume_fraction
      e. tortuosity_factor
      f. pores_characteristic_dimension
      g. pores_geometry_factor
      h. mass_density
      i. electrical_resistivity
      j. heat_capacity
      k. thermal_conductivity
    * solution_phase_X
      a. electrical_resistivity
      b. mass_density
    * metal_foil_X
      a. mass_density
      b. electrical_resistivity
      c. heat_capacity
      d. thermal_conductivity
    * inhomogeneous (bool)
    * parameters (unsigned int)
    * parameter_X (X in [0, parameters))
      a. path (matrix_phase_X/solution_phase_X/metal_foil_X.property)
      b. distribution_type (string)
      c. range (string)
      d. mean (double)
      e. standard_deviation (double)
      f. location (double)
      g. scale (double)
  6. solver
    * max_iter (unsigned int)
    * rel_tolerance (double)
    * abs_tolerance (double)
    * n_threads (unsigned int)

