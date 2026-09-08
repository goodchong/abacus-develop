#include "H_Hartree_pw.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "pot_local.h"
#include "pot_xc.h"
#include "potential_new.h"

namespace elecstate
{

PotBase* Potential::get_pot_type(const std::string& pot_type)
{
    ModuleBase::TITLE("Potential", "get_pot_type");
    if (pot_type == "local")
    {
        return new PotLocal(this->vloc_, &(this->structure_factor_->strucFac), this->rho_basis_, this->vl_of_0);
    }
    else if (pot_type == "hartree")
    {
        return new PotHartree(this->rho_basis_);
    }
    else if (pot_type == "xc")
    {
        return new PotXC(this->rho_basis_, this->etxc_, this->vtxc_, &(this->vofk_eff));
    }
    else
    {
        ModuleBase::WARNING_QUIT("Potential::get_pot_type", "Please input correct component of potential!");
        __builtin_unreachable();
    }
}

} // namespace elecstate
