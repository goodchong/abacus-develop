if (ENABLE_FLOAT_FFTW)
  list (APPEND FFT_SRC
    module_fft/fft_cpu_float.cpp
  )
endif()
if (USE_CUDA)
  list (APPEND FFT_SRC
    module_fft/fft_cuda.cpp
  )
endif()
if (USE_ROCM)
  list (APPEND FFT_SRC
    module_fft/fft_rocm.cpp
  )
endif()

list(APPEND objects
    pw_basis.cpp
    pw_basis_k.cpp
    pw_basis_sup.cpp
    pw_distributeg.cpp
    pw_distributeg_method1.cpp
    pw_distributeg_method2.cpp
    pw_distributer.cpp
    pw_init.cpp
    pw_transform.cpp
    pw_transform_k.cpp
    module_fft/fft_bundle.cpp
    module_fft/fft_cpu.cpp
    ${FFT_SRC}
)

add_library(
    planewave
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(planewave)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
    add_subdirectory(test_serial)
    add_subdirectory(kernels/test)
  endif()
endif()
