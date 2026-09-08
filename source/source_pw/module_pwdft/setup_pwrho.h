#ifndef SETUP_PWRHO_H
#define SETUP_PWRHO_H

#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

struct Input_para;

namespace pw
{

void setup_pwrho(UnitCell& ucell,
                 bool& allocated,
                 ModulePW::PW_Basis*& pw_rho,
                 ModulePW::PW_Basis*& pw_rhod,
                 ModulePW::PW_Basis_Big*& pw_big,
                 const Input_para& input);
void teardown_pwrho(bool& allocated,
                    ModulePW::PW_Basis*& pw_rho,
                    ModulePW::PW_Basis*& pw_rhod);

} // namespace pw

#endif
