set(BASE_CORE_CPU_SRCS
    base/core/refcount.cpp
    base/core/cpu_allocator.cpp
)

set(ATen_CPU_SRCS ${ATen_CPU_SRCS} ${BASE_CORE_CPU_SRCS} PARENT_SCOPE)
