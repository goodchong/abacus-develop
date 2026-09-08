

list(APPEND objects
    setup_pwrho.cpp
    vl_pw.cpp
    structure_factor.cpp
    soc.cpp
    parallel_grid.cpp
)

add_library(
    module_pwdft
    OBJECT
    ${objects}
)
