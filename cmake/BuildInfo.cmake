function(setup_build_info)
  if(NOT CMAKE_BUILD_TYPE)
    set(ABACUS_BUILD_TYPE "Custom")
  else()
    set(ABACUS_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
  endif()

  if(ENABLE_MPI)
    set(ABACUS_MPI_STATUS "yes")
  else()
    set(ABACUS_MPI_STATUS "no")
  endif()
  if(ENABLE_OPENMP)
    set(ABACUS_OPENMP_STATUS "yes")
  else()
    set(ABACUS_OPENMP_STATUS "no")
  endif()

  configure_file(
    "${CMAKE_SOURCE_DIR}/source/source_io/build_info.h.in"
    "${CMAKE_BINARY_DIR}/source/source_io/build_info.h"
    @ONLY)
endfunction()
