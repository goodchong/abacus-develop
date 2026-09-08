list(APPEND objects
    pw_basis.cpp
    pw_basis_sup.cpp
    pw_distributeg.cpp
    pw_distributeg_method1.cpp
    pw_distributeg_method2.cpp
    pw_distributer.cpp
    pw_init.cpp
    pw_transform.cpp
    ../../source_base/module_fft/fft_bundle.cpp
    ../../source_base/module_fft/fft_cpu.cpp
)

add_library(
    planewave
    OBJECT
    ${objects}
)
