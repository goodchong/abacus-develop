#ifndef LCAO_H0_BASIS_H
#define LCAO_H0_BASIS_H

#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/unitcell.h"

namespace LCAO_domain
{

void init_basis_lcao(Parallel_Orbitals& pv,
                     const double& lcao_ecut,
                     const double& lcao_dk,
                     UnitCell& ucell,
                     TwoCenterBundle& two_center_bundle,
                     LCAO_Orbitals& orb);

} // namespace LCAO_domain

#endif
