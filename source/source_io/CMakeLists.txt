# Note: Help data is now embedded in Input_Item objects in read_input_item_*.cpp
# The old input_help_data.h generation has been removed.

list(APPEND objects
    module_parameter/input_conv.cpp
    module_output/read_cube.cpp
    module_chgpot/rhog_io.cpp
    module_output/binstream.cpp
    module_output/sparse_matrix.cpp
    module_hs/write_h0.cpp
    parse_args.cpp
    input_help.cpp
    module_output/ucell_io.cpp
)

add_library(
    io_input
    OBJECT
    module_parameter/read_input_item_h0.cpp
    module_parameter/read_input.cpp
    module_parameter/read_set_globalv.cpp
)

add_library(
    io_basic
    OBJECT
    ${objects}
)

# Include directory for build_info.h (generated at build time)
target_include_directories(
    io_basic
    PUBLIC
    ${CMAKE_BINARY_DIR}/source/source_io
)
