#include "lcao_h0_basis.h"

#include "source_base/global_function.h"
#include "source_base/parallel_comm.h"
#include "source_base/tool_title.h"
#include "source_io/module_parameter/parameter.h"
#include "LCAO_nonlocal_info.h"

namespace LCAO_domain
{

void init_basis_lcao(Parallel_Orbitals& pv,
        const double &lcao_ecut,
        const double &lcao_dk,
		UnitCell& ucell,
        TwoCenterBundle& two_center_bundle,
        LCAO_Orbitals& orb
)
{
    ModuleBase::TITLE("H0Runner", "init_basis_lcao");

    const int nlocal = PARAM.globalv.nlocal;
    int nb2d = PARAM.inp.nb2d;
    // autoset NB2D first
    if (nb2d == 0)
    {
        if (nlocal > 0)
        {
            nb2d = (PARAM.inp.nspin == 4) ? 2 : 1;
        }
        if (nlocal > 500)
        {
            nb2d = 32;
        }
        if (nlocal > 1000)
        {
            nb2d = 64;
        }
    }

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.

    two_center_bundle.build_orb(ucell.ntype, ucell.orbital_fn.data(), PARAM.inp.orbital_dir);
    two_center_bundle.to_LCAO_Orbitals(orb, lcao_ecut, lcao_dk, false);

    auto* lcao_nl = new LCAONonlocalInfo();
    lcao_nl->setupNonlocal(ucell.ntype,
                           ucell.atoms,
                           GlobalV::ofs_running,
                           orb,
                           PARAM.inp.lspinorb,
                           PARAM.inp.nspin);
    ucell.infoNL.reset(lcao_nl);
    two_center_bundle.build_beta(ucell.ntype, lcao_nl->get_nonlocal().Beta);

    two_center_bundle.tabulate();

    // setup_2d_division
#ifdef __MPI
    // storage form of H and S matrices on each processor
    // is determined in 'divide_HS_2d' subroutine

    int try_nb = pv.init(nlocal, nlocal, nb2d, DIAG_WORLD);
    if (try_nb != 0)
    {
        // fall back to the minimum size, 1 or 2 (nspin=4)
        const int min_size = (PARAM.inp.nspin == 4) ? 2 : 1;
        pv.set(nlocal, nlocal, min_size, pv.blacs_ctxt);
    }

#else
    pv.set_serial(nlocal, nlocal);
#endif

    pv.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, nlocal);

    return;
}

}
