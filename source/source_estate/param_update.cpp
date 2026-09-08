#include "param_update.h"
#include "source_io/module_parameter/parameter.h"

namespace elecstate {

void ParamUpdater::update_from_atoms_info(const AtomsInfoResult& atoms_info)
{
    PARAM.input.nelec = atoms_info.nelec;
    PARAM.input.nupdown = atoms_info.nupdown;
    PARAM.sys.nlocal = atoms_info.nlocal;
}

}
