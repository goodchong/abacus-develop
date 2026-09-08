list(APPEND objects
    module_pot/H_Hartree_pw.cpp
    module_pot/pot_xc.cpp
    module_pot/pot_local.cpp
    module_pot/potential_new.cpp
    module_pot/potential_types.cpp
    module_charge/charge.cpp
    module_charge/charge_init.cpp
    cal_ux.cpp
    read_pseudo.cpp
    param_update.cpp
    basis_index.cpp
)

add_library(
    elecstate
    OBJECT
    ${objects}
)
