find_package(getkw CONFIG QUIET)
if(TARGET getkw::getkw)
  get_property(_loc TARGET getkw::getkw PROPERTY LOCATION)
  message(STATUS "Found getkw: ${_loc} (found version ${getkw_VERSION})")
else()
  message(STATUS "Suitable getkw could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Populate(getkw_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/dev-cafe/libgetkw.git
    GIT_TAG
      bb3eb9726bc443d34f89629b025e90a4fc793c73 # Preferable to have a tag for a release
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DENABLE_TESTS=TRUE
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    )

  add_subdirectory(
    ${getkw_sources_SOURCE_DIR}
    ${getkw_sources_BINARY_DIR}
    )
endif()