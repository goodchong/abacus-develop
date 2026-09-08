#include "source_basis/module_nao/two_center_bundle.h"

#include "source_base/global_variable.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_common.h"
#include "source_base/ylm.h"
#include "source_basis/module_nao/real_gaunt_table.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

void TwoCenterBundle::build_orb(int ntype,
                                const std::string* file_orb0,
                                const std::string& orbital_dir)
{
    std::vector<std::string> file_orb(ntype);
    if (GlobalV::MY_RANK == 0)
    {
        std::transform(file_orb0,
                       file_orb0 + ntype,
                       file_orb.begin(),
                       [&orbital_dir](const std::string& file) {
                           return orbital_dir + file;
                       });
    }
#ifdef __MPI
    Parallel_Common::bcast_string(file_orb.data(), ntype);
#endif

    orb_.reset(new RadialCollection);
    orb_->build(ntype, file_orb.data());
}

void TwoCenterBundle::build_beta(int ntype, Numerical_Nonlocal* nl)
{
    beta_.reset(new RadialCollection);
    beta_->build(ntype, nl);
}

void TwoCenterBundle::tabulate()
{
    ModuleBase::SphericalBesselTransformer sbt(true);
    orb_->set_transformer(sbt);
    if (beta_)
    {
        beta_->set_transformer(sbt);
    }

    double rmax = orb_->rcut_max();
    if (beta_)
    {
        rmax = std::max(rmax, beta_->rcut_max());
    }
    const double dr = 0.01;
    const double cutoff = 2.0 * rmax;
    const int nr = static_cast<int>(rmax / dr) + 1;

    orb_->set_uniform_grid(true, nr, cutoff, 'i', true);
    if (beta_)
    {
        beta_->set_uniform_grid(true, nr, cutoff, 'i', true);
    }

    kinetic_orb.reset(new TwoCenterIntegrator);
    kinetic_orb->tabulate(*orb_, *orb_, 'T', nr, cutoff);
    ModuleBase::Memory::record("TwoCenterTable: Kinetic",
                               kinetic_orb->table_memory());

    if (beta_)
    {
        overlap_orb_beta.reset(new TwoCenterIntegrator);
        overlap_orb_beta->tabulate(*orb_, *beta_, 'S', nr, cutoff);
        ModuleBase::Memory::record("TwoCenterTable: Nonlocal",
                                   overlap_orb_beta->table_memory());
    }

    ModuleBase::Memory::record("RealGauntTable",
                               RealGauntTable::instance().memory());
    sbt.clear();
}

void TwoCenterBundle::to_LCAO_Orbitals(LCAO_Orbitals& orbitals,
                                       const double lcao_ecut,
                                       const double lcao_dk,
                                       const bool out_element_info) const
{
    orbitals.ntype = orb_->ntype();
    orbitals.rcutmax_Phi = orb_->rcut_max();
    orbitals.dr_uniform = 0.001;
    orbitals.dk = lcao_dk;

    if (lcao_ecut < 20)
    {
        orbitals.kmesh
            = static_cast<int>(2 * std::sqrt(lcao_ecut) / orbitals.dk) + 4;
    }
    else
    {
        orbitals.kmesh
            = static_cast<int>(std::sqrt(lcao_ecut) / orbitals.dk) + 4;
    }
    orbitals.kmesh += 1 - orbitals.kmesh % 2;

    delete[] orbitals.Phi;
    orbitals.Phi = new Numerical_Orbital[orb_->ntype()];
    for (int itype = 0; itype < orb_->ntype(); ++itype)
    {
        (*orb_)(itype).to_numerical_orbital(orbitals.Phi[itype],
                                             orbitals.kmesh,
                                             orbitals.dk,
                                             out_element_info,
                                             false);
    }
}
