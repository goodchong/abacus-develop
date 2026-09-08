#ifndef POTENTIAL_NEW_H
#define POTENTIAL_NEW_H

#include "pot_base.h"
#include "source_pw/module_pwdft/structure_factor.h"

#include <vector>

namespace elecstate
{

class Potential : public PotBase
{
  public:
    Potential(const ModulePW::PW_Basis* rho_basis,
              const UnitCell* ucell,
              const ModuleBase::matrix* vloc,
              Structure_Factor* structure_factor,
              double* etxc,
              double* vtxc);
    ~Potential();

    void pot_register(const std::vector<std::string>& component_names);
    void init_pot(const Charge* charge);
    PotBase* get_pot_type(const std::string& component_name);

    const double* get_eff_v(int spin) const
    {
        return this->v_eff.nc > 0 ? &this->v_eff(spin, 0) : nullptr;
    }

  private:
    void allocate();
    void update_from_charge(const Charge* charge);
    void cal_v_eff(const Charge* charge,
                   const UnitCell* ucell,
                   ModuleBase::matrix& v_eff) override;
    void cal_fixed_v(double* local_potential) override;

    std::vector<double> v_eff_fixed;
    ModuleBase::matrix v_eff;
    ModuleBase::matrix vofk_eff;
    bool fixed_done = false;
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;
    double vl_of_0 = 0.0;
    std::vector<PotBase*> components;
    const UnitCell* ucell_ = nullptr;
    const ModuleBase::matrix* vloc_ = nullptr;
    Structure_Factor* structure_factor_ = nullptr;
};

} // namespace elecstate

#endif
