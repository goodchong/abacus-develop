add_subdirectory(module_hcontainer)

list(APPEND objects
        module_operator_lcao/ekinetic.cpp
        module_operator_lcao/nonlocal.cpp
        LCAO_init_basis.cpp
        setup_nonlocal.cpp
        record_adj.cpp
    )

add_library(
    hamilt_lcao
    OBJECT
    ${objects})
