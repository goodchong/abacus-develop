#include "pot_xc.h"

#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

namespace elecstate
{

void PotXC::cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff)
{
    ModuleBase::TITLE("PotXC", "cal_veff");
    ModuleBase::timer::start("PotXC", "cal_veff");
    const int nrxx_current = chg->nrxx;
    
    //----------------------------------------------------------
    //  calculate the exchange-correlation potential
    //----------------------------------------------------------

    const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
        = XC_Functional::v_xc(nrxx_current,
                              chg,
                              ucell,
                              PARAM.inp.nspin,
                              PARAM.globalv.domag,
                              PARAM.globalv.domag_z);
    *(this->etxc_) = std::get<0>(etxc_vtxc_v);
    *(this->vtxc_) = std::get<1>(etxc_vtxc_v);
    v_eff += std::get<2>(etxc_vtxc_v);
    ModuleBase::timer::end("PotXC", "cal_veff");
}

} // namespace elecstate
