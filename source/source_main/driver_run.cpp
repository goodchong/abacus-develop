#include "source_main/driver.h"

#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_cell/check_atomic_stru.h"
#include "source_h0/h0_runner.h"
#include "source_io/module_parameter/parameter.h"

void Driver::driver_run()
{
    ModuleBase::TITLE("Driver", "driver_run");

#ifndef __LCAO
    ModuleBase::WARNING_QUIT("Driver", "The H0-only executable must be built with LCAO support.");
#endif

    UnitCell ucell;
    ucell.setup(PARAM.inp.ntype);
    ucell.setup_cell(PARAM.globalv.global_in_stru,
                     GlobalV::ofs_running,
                     PARAM.inp.nspin,
                     PARAM.inp.orbital_dir,
                     PARAM.inp.noncolin);
    unitcell::check_atomic_stru(ucell, PARAM.inp.min_dist_coef);

    ModuleH0::H0Runner runner;
    runner.initialize(ucell, PARAM.inp);
    runner.run(ucell);
    runner.finalize();
}
