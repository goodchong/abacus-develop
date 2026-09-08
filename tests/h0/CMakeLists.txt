find_package(Python3 COMPONENTS Interpreter REQUIRED)

foreach(_case IN ITEMS si2 fe2_collinear si_soc)
  foreach(_h0_type IN ITEMS core full)
    add_test(
      NAME "H0_${_case}_${_h0_type}"
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
              --abacus "$<TARGET_FILE:abacus>"
              --source-root "${PROJECT_SOURCE_DIR}"
              --case "${_case}"
              --h0-type "${_h0_type}")
    set_tests_properties("H0_${_case}_${_h0_type}" PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1"
      LABELS "H0")
  endforeach()
endforeach()

add_test(
  NAME H0_init_chg_auto_fallback
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type full --init-chg auto)

add_test(
  NAME H0_init_chg_file
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type full --init-chg file
          --density-file "${CMAKE_CURRENT_SOURCE_DIR}/reference/si2_chg.cube"
          --skip-baseline)

add_test(
  NAME H0_init_chg_file_missing
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type full --init-chg file
          --expect-failure "Failed to read in charge density from file")

add_test(
  NAME H0_reject_legacy_kpoint_switch
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type core --extra-input "gamma_only 1"
          --expect-failure "THE PARAMETER NAME 'gamma_only' IS INCORRECT")

add_test(
  NAME H0_require_h0_type
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type core --remove-input h0_type
          --expect-failure "h0_type must be explicitly set to core or full")

add_test(
  NAME H0_reject_scf_calculation
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type core --set-input calculation=scf
          --expect-failure "calculation must be get_h0")

add_test(
  NAME H0_reject_pw_basis
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type core --set-input basis_type=pw
          --expect-failure "basis_type must be lcao")

add_test(
  NAME H0_reject_nspin_three
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type core --set-input nspin=3
          --expect-failure "nspin must be 1, 2, or 4")

add_test(
  NAME H0_reject_meta_gga
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type full --set-input dft_functional=scan
          --expect-failure "dft_functional must select a built-in LDA or GGA functional")

add_test(
  NAME H0_threshold_strict_boundary
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_threshold.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}")

add_test(
  NAME H0_precision_lower_bound
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si_soc --h0-type core --skip-baseline
          --extra-input "h0_precision 1")

add_test(
  NAME H0_precision_upper_bound_rejected
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si_soc --h0-type core --extra-input "h0_precision 18"
          --expect-failure "h0_precision must be between 1 and 17")

add_test(
  NAME H0_openmp_four_threads
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
          --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
          --case si2 --h0-type full --omp-threads 4)

set_tests_properties(
  H0_init_chg_auto_fallback
  H0_init_chg_file
  H0_init_chg_file_missing
  H0_reject_legacy_kpoint_switch
  H0_require_h0_type
  H0_reject_scf_calculation
  H0_reject_pw_basis
  H0_reject_nspin_three
  H0_reject_meta_gga
  H0_threshold_strict_boundary
  H0_precision_lower_bound
  H0_precision_upper_bound_rejected
  H0_openmp_four_threads
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "H0")

if(ENABLE_MPI)
  foreach(_ranks IN ITEMS 1 2 4)
    foreach(_case IN ITEMS si2 fe2_collinear si_soc)
      foreach(_h0_type IN ITEMS core full)
        add_test(
          NAME "H0_mpi${_ranks}_${_case}_${_h0_type}"
          COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_h0.py"
                  --abacus "$<TARGET_FILE:abacus>" --source-root "${PROJECT_SOURCE_DIR}"
                  --case "${_case}" --h0-type "${_h0_type}"
                  --mpiexec "${MPIEXEC_EXECUTABLE}" --mpi-ranks "${_ranks}")
        set_tests_properties("H0_mpi${_ranks}_${_case}_${_h0_type}" PROPERTIES
          ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "H0;MPI" PROCESSORS "${_ranks}")
      endforeach()
    endforeach()
  endforeach()
endif()
