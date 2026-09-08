list(APPEND objects
    base_matrix.cpp
    atom_pair.cpp
    hcontainer.cpp
    output_hcontainer.cpp
    transfer.cpp
    func_transfer.cpp
)

add_library(
    hcontainer
    OBJECT
    ${objects}
)
