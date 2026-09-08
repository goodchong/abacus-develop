#include "electron_count.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"

namespace unitcell
{

void set_default_electron_count(const Atom* atoms,
                                const int ntype,
                                double& nelec)
{
    ModuleBase::TITLE("UnitCell", "set_default_electron_count");
    if (nelec != 0.0)
    {
        return;
    }
    for (int it = 0; it < ntype; ++it)
    {
        const double species_electrons = atoms[it].ncpp.zv * atoms[it].na;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                    "Electron number of element "
                                        + atoms[it].label,
                                    atoms[it].ncpp.zv);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                    "Total electron number of element "
                                        + atoms[it].label,
                                    species_electrons);
        nelec += species_electrons;
    }
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "Autoset the number of electrons",
                                nelec);
}

} // namespace unitcell
