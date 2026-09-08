#include "source_io/module_output/ucell_io.h"
#include "source_base/constants.h"

namespace ModuleIO {

void UcellIO::write_ucell(std::ofstream& ofs, const UnitCell* ucell)
{
    // write the UnitCell information
    ofs << " " << ucell->latName << std::endl;
    ofs << " " << ucell->lat0 * ModuleBase::BOHR_TO_A << std::endl;
    ofs << " " << ucell->latvec.e11 << " " << ucell->latvec.e12 << " " << ucell->latvec.e13 << std::endl;
    ofs << " " << ucell->latvec.e21 << " " << ucell->latvec.e22 << " " << ucell->latvec.e23 << std::endl;
    ofs << " " << ucell->latvec.e31 << " " << ucell->latvec.e32 << " " << ucell->latvec.e33 << std::endl;
    for (int it = 0; it < ucell->ntype; it++)
    {
        ofs << " " << ucell->atoms[it].label;
    }
    ofs << std::endl;
    for (int it = 0; it < ucell->ntype; it++)
    {
        ofs << " " << ucell->atoms[it].na;
    }
    ofs << std::endl;
    ofs << " Direct" << std::endl;
    for (int it = 0; it < ucell->ntype; it++)
    {
        Atom* atom = &ucell->atoms[it];
        ofs << std::setprecision(15);
        for (int ia = 0; ia < ucell->atoms[it].na; ia++)
        {
            ofs << " " << atom->taud[ia].x << " " << atom->taud[ia].y << " " << atom->taud[ia].z << std::endl;
        }
    }
}

} // namespace ModuleIO
