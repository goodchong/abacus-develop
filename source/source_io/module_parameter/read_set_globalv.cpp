#include "read_input.h"
#include "read_input_tool.h"

#include "source_base/parallel_common.h"

#include <cmath>

namespace ModuleIO
{

void ReadInput::set_globalv(const Input_para& inp, System_para& sys)
{
    sys.npol = inp.nspin == 4 ? 2 : 1;
    sys.domag = inp.nspin == 4 && inp.noncolin;
    sys.domag_z = inp.nspin == 4 && !inp.noncolin;
    sys.two_fermi = inp.nspin == 2 && inp.nupdown != 0.0;
}

void ReadInput::set_global_dir(const Input_para& inp, System_para& sys)
{
    sys.global_out_dir = to_dir("OUT." + inp.suffix);
    sys.global_readin_dir = to_dir(inp.read_file_dir);
    sys.global_in_stru = inp.stru_file;
    sys.log_file = "running_get_h0.log";

#ifdef __MPI
    Parallel_Common::bcast_string(sys.global_in_card);
    Parallel_Common::bcast_string(sys.global_out_dir);
    Parallel_Common::bcast_string(sys.global_readin_dir);
    Parallel_Common::bcast_string(sys.global_in_stru);
#endif
}

} // namespace ModuleIO
