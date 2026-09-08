#include "basis_index.h"

#include "source_base/tool_title.h"

#include <algorithm>
#include <cassert>

namespace elecstate
{

void setup_lcao_basis_indices(UnitCell& ucell,
                              Atom* atoms,
                              const int nspin,
                              const int nlocal,
                              const int npol)
{
    ModuleBase::TITLE("UnitCell", "setup_lcao_basis_indices");
    assert(ucell.ntype > 0);
    assert(ucell.nat > 0);

    ucell.namax = 0;
    int local_basis_count = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        atoms[it].set_index();
        ucell.namax = std::max(atoms[it].na, ucell.namax);
        const int species_basis = atoms[it].nw * atoms[it].na;
        local_basis_count += nspin == 4 ? 2 * species_basis : species_basis;
    }
    assert(local_basis_count == nlocal);

    ucell.itia2iat.create(ucell.ntype, ucell.namax);
    ucell.set_iat2iwt(npol);
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < atoms[it].na; ++ia)
        {
            ucell.itia2iat(it, ia) = iat++;
        }
    }
}

void set_maximum_pseudo_mesh(int& meshx,
                             const Atom* atoms,
                             const int ntype)
{
    meshx = 0;
    for (int it = 0; it < ntype; ++it)
    {
        meshx = std::max(meshx, atoms[it].ncpp.msh);
    }
}

} // namespace elecstate
