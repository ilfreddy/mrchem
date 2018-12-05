set(_h2_pol_lda_launcher " ")

if(ENABLE_MPI)
    set(_h2_pol_lda_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME "H2_static_polarizability_LDA"
  LABELS "mrchem;h2_pol_lda"
  COST 100
  LAUNCH_AGENT ${_h2_pol_lda_launcher}
  )
