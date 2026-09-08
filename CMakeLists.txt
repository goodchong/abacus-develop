cmake_minimum_required(VERSION 3.16)

project(
  ABACUS_H0
  VERSION 3.11.0
  DESCRIPTION "One-shot LCAO initial-Hamiltonian matrix generator"
  LANGUAGES CXX)

option(ENABLE_MPI "Enable MPI and ScaLAPACK distributed H(R) construction" ON)
option(ENABLE_OPENMP "Enable OpenMP" ON)
option(BUILD_TESTING "Build the H0 regression tests" OFF)
option(COMMIT_INFO "Embed the current Git revision" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
set(ABACUS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(ABACUS_BIN_NAME abacus)

set_property(GLOBAL PROPERTY ABACUS_FEATURE_DEFINITIONS "")

function(abacus_normalize_definitions out_var)
  set(_defs)
  foreach(_def IN LISTS ARGN)
    if(_def MATCHES "^-D(.+)")
      list(APPEND _defs "${CMAKE_MATCH_1}")
    else()
      list(APPEND _defs "${_def}")
    endif()
  endforeach()
  set(${out_var} ${_defs} PARENT_SCOPE)
endfunction()

function(abacus_add_feature_definitions)
  abacus_normalize_definitions(_defs ${ARGN})
  set_property(GLOBAL APPEND PROPERTY ABACUS_FEATURE_DEFINITIONS ${_defs})
endfunction()

abacus_add_feature_definitions(__LCAO __FFTW3 USE_NEW_TWO_CENTER)

if(COMMIT_INFO)
  find_package(Git QUIET)
  if(Git_FOUND)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:%h
      OUTPUT_VARIABLE GIT_COMMIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      RESULT_VARIABLE GIT_HASH_RESULT)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" log -1 --format=%cd
      OUTPUT_VARIABLE GIT_COMMIT_DATE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      RESULT_VARIABLE GIT_DATE_RESULT)
    if(GIT_HASH_RESULT EQUAL 0 AND GIT_DATE_RESULT EQUAL 0)
      file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/commit.h"
           "#define COMMIT \"${GIT_COMMIT_HASH} (${GIT_COMMIT_DATE})\"\n")
      set(ABACUS_COMMIT_INFO_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
      abacus_add_feature_definitions(COMMIT_INFO)
    endif()
  endif()
endif()

find_package(Threads REQUIRED)
find_package(Lapack REQUIRED)
find_package(FFTW3 REQUIRED)

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if(ENABLE_MPI)
  find_package(MPI COMPONENTS CXX REQUIRED)
  find_package(ScaLAPACK REQUIRED)
  abacus_add_feature_definitions(__MPI)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildInfo.cmake")
setup_build_info()

add_subdirectory(source)

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests/h0)
endif()

install(TARGETS abacus RUNTIME DESTINATION bin)
