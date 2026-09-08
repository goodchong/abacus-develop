add_library(abacus_compile_requirements INTERFACE)
add_library(abacus::compile_requirements ALIAS abacus_compile_requirements)
target_include_directories(abacus_compile_requirements INTERFACE
  "${ABACUS_SOURCE_DIR}"
  "${ABACUS_SOURCE_DIR}/source_base/module_container")
if(ABACUS_COMMIT_INFO_INCLUDE_DIR)
  target_include_directories(abacus_compile_requirements INTERFACE
    "${ABACUS_COMMIT_INFO_INCLUDE_DIR}")
endif()

get_property(_h0_definitions GLOBAL PROPERTY ABACUS_FEATURE_DEFINITIONS)
target_compile_definitions(abacus_compile_requirements INTERFACE ${_h0_definitions})

if(ENABLE_MPI)
  target_link_libraries(abacus_compile_requirements INTERFACE MPI::MPI_CXX)
endif()
if(ENABLE_OPENMP)
  target_link_libraries(abacus_compile_requirements INTERFACE OpenMP::OpenMP_CXX)
endif()

add_executable(abacus source_main/main.cpp)
set_target_properties(abacus PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_subdirectory(source_base)
add_subdirectory(source_cell)
add_subdirectory(source_estate)
add_subdirectory(source_hamilt/module_xc)
add_subdirectory(source_pw)
add_subdirectory(source_lcao)
add_subdirectory(source_basis/module_ao)
add_subdirectory(source_basis/module_nao)
add_subdirectory(source_basis/module_pw)
add_subdirectory(source_h0)
add_subdirectory(source_lcao/module_gint)
add_subdirectory(source_io)
add_subdirectory(source_io/module_parameter)

add_library(driver OBJECT
  source_main/driver.cpp
  source_main/driver_run.cpp)

# CPU implementations required by the retained matrix, FFT and tensor support.
add_library(device OBJECT
  source_base/module_device/device.cpp
  source_base/module_device/device_helpers.cpp
  source_base/module_device/output_device.cpp
  source_base/module_device/memory_op.cpp
  source_base/kernels/math_kernel_op.cpp
  source_base/kernels/math_kernel_op_vec.cpp
  source_base/kernels/math_ylm_op.cpp)
target_link_libraries(base PUBLIC device)

set(_h0_numerical_libraries
  FFTW3::FFTW3
  LAPACK::LAPACK
  BLAS::BLAS
  Threads::Threads)
if(ENABLE_OPENMP)
  list(PREPEND _h0_numerical_libraries FFTW3::FFTW3_OMP OpenMP::OpenMP_CXX)
endif()
if(ENABLE_MPI)
  list(PREPEND _h0_numerical_libraries ScaLAPACK::ScaLAPACK MPI::MPI_CXX)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  list(APPEND _h0_numerical_libraries gfortran)
endif()
if(NOT MSVC)
  list(APPEND _h0_numerical_libraries m)
endif()

target_link_libraries(abacus PRIVATE
  driver
  h0_runner
  elecstate
  module_pwdft
  xc_
  io_basic
  io_input
  neighbor
  planewave
  cell
  parameter
  container
  device
  base
  hamilt_lcao
  orb
  gint
  hcontainer
  numerical_atomic_orbitals
  ${_h0_numerical_libraries})

function(h0_apply_requirements directory)
  get_property(_targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_target IN LISTS _targets)
    get_target_property(_type "${_target}" TYPE)
    if(NOT _type STREQUAL "INTERFACE_LIBRARY" AND NOT _type STREQUAL "UTILITY")
      target_link_libraries("${_target}" PRIVATE abacus::compile_requirements)
    endif()
  endforeach()
  get_property(_subdirs DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
  foreach(_subdir IN LISTS _subdirs)
    string(FIND "${_subdir}/" "${PROJECT_SOURCE_DIR}/source/" _inside_source)
    if(_inside_source EQUAL 0)
      h0_apply_requirements("${_subdir}")
    endif()
  endforeach()
endfunction()

h0_apply_requirements("${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(abacus PRIVATE abacus::compile_requirements)
