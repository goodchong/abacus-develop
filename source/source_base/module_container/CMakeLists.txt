# FetchContent_MakeAvailable was not added until CMake 3.14
cmake_minimum_required(VERSION 3.14)

set(ATen_CPU_SRCS)
set(ATen_CPU_DEPENDENCY_LIBS)

add_subdirectory(base/core)

add_subdirectory(ATen/core)

# The H0-only program uses Tensor as CPU-owned storage for radial tables.  It
# needs allocation/copy kernels, but no tensor algebra or LAPACK eigensolvers.
list(APPEND ATen_CPU_SRCS ATen/kernels/memory_impl.cpp)

add_library(container STATIC ${ATen_CPU_SRCS})

target_link_libraries(container PUBLIC ${ATen_CPU_DEPENDENCY_LIBS})
