set(_li_pol_lda_launcher " ")

if(ENABLE_MPI)
    set(_li_pol_lda_launcher "mpirun -np 2")
endif()

add_integration_test(
  NAME "Li_static_polarizability_LDA"
  LABELS "mrchem;li_pol_lda"
  COST 100
  LAUNCH_AGENT ${_li_pol_lda_launcher}
  )
