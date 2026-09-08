#include "input_conv.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

namespace Input_Conv
{

void Convert()
{
    ModuleBase::TITLE("Input_Conv", "Convert");
    ModuleBase::timer::start("Input_Conv", "Convert");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pseudo_dir", PARAM.inp.pseudo_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "orbital_dir", PARAM.inp.orbital_dir);

    // The H0 workflow has no k-point pools, but the retained FFT and LCAO
    // distribution code expects a single initialized pool.
    GlobalV::KPAR = 1;
    ModuleBase::timer::end("Input_Conv", "Convert");
}

} // namespace Input_Conv
