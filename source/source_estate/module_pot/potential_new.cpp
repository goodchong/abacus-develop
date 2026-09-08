#include "potential_new.h"

#include "source_base/global_function.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

namespace elecstate
{

Potential::Potential(const ModulePW::PW_Basis* rho_basis,
                     const UnitCell* ucell,
                     const ModuleBase::matrix* vloc,
                     Structure_Factor* structure_factor,
                     double* etxc,
                     double* vtxc)
    : etxc_(etxc),
      vtxc_(vtxc),
      ucell_(ucell),
      vloc_(vloc),
      structure_factor_(structure_factor)
{
    this->rho_basis_ = rho_basis;
    this->fixed_mode = true;
    this->dynamic_mode = true;
    this->allocate();
}

Potential::~Potential()
{
    for (std::size_t i = 0; i < this->components.size(); ++i)
    {
        delete this->components[i];
    }
}

void Potential::pot_register(const std::vector<std::string>& component_names)
{
    for (std::size_t i = 0; i < this->components.size(); ++i)
    {
        delete this->components[i];
    }
    this->components.clear();
    for (std::size_t i = 0; i < component_names.size(); ++i)
    {
        this->components.push_back(this->get_pot_type(component_names[i]));
    }
    this->fixed_done = false;
}

void Potential::allocate()
{
    const int nspin = PARAM.inp.nspin;
    const int nrxx = this->rho_basis_->nrxx;
    this->v_eff_fixed.resize(nrxx);
    this->v_eff.create(nspin, nrxx);
    ModuleBase::Memory::record("Pot::veff_fix", sizeof(double) * nrxx);
    ModuleBase::Memory::record("Pot::veff", sizeof(double) * nspin * nrxx);
}

void Potential::init_pot(const Charge* charge)
{
    this->fixed_done = false;
    this->update_from_charge(charge);
}

void Potential::update_from_charge(const Charge* charge)
{
    if (!this->fixed_done)
    {
        this->cal_fixed_v(this->v_eff_fixed.data());
        this->fixed_done = true;
    }
    this->cal_v_eff(charge, this->ucell_, this->v_eff);
}

void Potential::cal_fixed_v(double* local_potential)
{
    ModuleBase::timer::start("Potential", "cal_fixed_v");
    this->v_eff_fixed.assign(this->v_eff_fixed.size(), 0.0);
    for (std::size_t i = 0; i < this->components.size(); ++i)
    {
        if (this->components[i]->fixed_mode)
        {
            this->components[i]->cal_fixed_v(local_potential);
        }
    }
    ModuleBase::timer::end("Potential", "cal_fixed_v");
}

void Potential::cal_v_eff(const Charge* charge,
                          const UnitCell* ucell,
                          ModuleBase::matrix& v_eff)
{
    ModuleBase::timer::start("Potential", "cal_veff");
    const int nspin = v_eff.nr;
    const int nrxx = v_eff.nc;
    v_eff.zero_out();
    for (int is = 0; is < nspin; ++is)
    {
        if (is == 0 || nspin == 2)
        {
            ModuleBase::GlobalFunc::COPYARRAY(this->v_eff_fixed.data(),
                                              &v_eff(is, 0),
                                              nrxx);
        }
    }
    for (std::size_t i = 0; i < this->components.size(); ++i)
    {
        if (this->components[i]->dynamic_mode)
        {
            this->components[i]->cal_v_eff(charge, ucell, v_eff);
        }
    }
    ModuleBase::timer::end("Potential", "cal_veff");
}

} // namespace elecstate
