#include "setup_pwrho.h"

#include "source_base/parallel_comm.h"
#include "source_io/module_parameter/input_parameter.h"

namespace pw
{

void setup_pwrho(UnitCell& ucell,
                 bool& allocated,
                 ModulePW::PW_Basis*& pw_rho,
                 ModulePW::PW_Basis*& pw_rhod,
                 ModulePW::PW_Basis_Big*& pw_big,
                 const Input_para& input)
{
    ModuleBase::TITLE("pw", "setup_pwrho");
    pw_rho = new ModulePW::PW_Basis_Big("cpu", "double");
    allocated = true;
    pw_rhod = pw_rho;
    pw_big = static_cast<ModulePW::PW_Basis_Big*>(pw_rho);
    pw_big->setbxyz(input.bx, input.by, input.bz);

#ifdef __MPI
    pw_rho->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    if (input.nx * input.ny * input.nz == 0)
    {
        pw_rho->initgrids(ucell.lat0, ucell.latvec, 4.0 * input.ecutwfc);
    }
    else
    {
        pw_rho->initgrids(ucell.lat0, ucell.latvec, input.nx, input.ny, input.nz);
    }
    pw_rho->initparameters(false, 4.0 * input.ecutwfc);
    pw_rho->fft_bundle.initfftmode(0);
    pw_rho->setuptransform();
    pw_rho->collect_local_pw();
    pw_rho->collect_uniqgg();

    GlobalV::ofs_running << " FFT grid: " << pw_rho->nx << " " << pw_rho->ny << " "
                         << pw_rho->nz << "; local real-space points: " << pw_rho->nrxx
                         << std::endl;
}

void teardown_pwrho(bool& allocated,
                    ModulePW::PW_Basis*& pw_rho,
                    ModulePW::PW_Basis*& pw_rhod)
{
    if (allocated)
    {
        delete pw_rho;
        pw_rho = nullptr;
        pw_rhod = nullptr;
        allocated = false;
    }
}

} // namespace pw
