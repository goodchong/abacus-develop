#ifndef CAL_ATOMS_INFO_H
#define CAL_ATOMS_INFO_H

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/electron_count.h"

struct AtomsInfoResult
{
    int nlocal = 0;
    double nelec = 0.0;
    double nupdown = 0.0;
};

class CalAtomsInfo
{
  public:
    AtomsInfoResult cal_atoms_info(Atom* atoms,
                                   int ntype,
                                   int nspin,
                                   bool two_fermi,
                                   double nelec,
                                   double nupdown) const
    {
        AtomsInfoResult result;
        if (nspin == 2 && !two_fermi)
        {
            for (int it = 0; it < ntype; ++it)
            {
                for (int ia = 0; ia < atoms[it].na; ++ia)
                {
                    result.nupdown += atoms[it].mag[ia];
                }
            }
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                        "The readin total magnetization",
                                        result.nupdown);
        }
        else if (nspin == 2)
        {
            result.nupdown = nupdown;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                        "The user-specified total magnetization",
                                        result.nupdown);
        }

        for (int it = 0; it < ntype; ++it)
        {
            atoms[it].set_index();
            const int species_basis = atoms[it].nw * atoms[it].na;
            result.nlocal += nspin == 4 ? 2 * species_basis : species_basis;
        }

        result.nelec = nelec;
        unitcell::set_default_electron_count(atoms, ntype, result.nelec);
        return result;
    }
};

#endif
