set(_h_el_field_launcher " ")

if(ENABLE_MPI)
    set(_h_el_field_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME "H_SCF_finite_electric_field"
  LABELS "mrchem;h_el_field"
  COST 100
  LAUNCH_AGENT ${_h_el_field_launcher}
  )
