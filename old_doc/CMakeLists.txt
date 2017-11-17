project(MRCPPdoc)
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

add_custom_target(doc 
  COMMAND doxygen
  WORKING_DIRECTORY ${MRCPP_BINARY_DIR}/doc
  )
