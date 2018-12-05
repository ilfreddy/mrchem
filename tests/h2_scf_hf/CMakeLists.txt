set(_h2_scf_hf_launcher " ")

if(ENABLE_MPI)
    set(_h2_scf_hf_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME "H2_SCF_HartreeFock"
  LABELS "mrchem;h2_scf_hf"
  COST 100
  LAUNCH_AGENT ${_h2_scf_hf_launcher}
  )
