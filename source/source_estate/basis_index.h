#ifndef BASIS_INDEX_H
#define BASIS_INDEX_H

#include "source_cell/unitcell.h"

namespace elecstate
{

void setup_lcao_basis_indices(UnitCell& ucell,
                              Atom* atoms,
                              int nspin,
                              int nlocal,
                              int npol);

void set_maximum_pseudo_mesh(int& meshx, const Atom* atoms, int ntype);

} // namespace elecstate

#endif
