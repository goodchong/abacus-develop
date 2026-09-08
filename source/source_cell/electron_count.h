#ifndef ELECTRON_COUNT_H
#define ELECTRON_COUNT_H

#include "source_cell/atom_spec.h"

namespace unitcell
{

void set_default_electron_count(const Atom* atoms, int ntype, double& nelec);

} // namespace unitcell

#endif
