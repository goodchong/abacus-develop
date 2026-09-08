#ifndef STRUCTURE_FACTOR_H
#define STRUCTURE_FACTOR_H

#include "source_base/complexmatrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/parallel_grid.h"

class Structure_Factor
{
  public:
    void set(const ModulePW::PW_Basis* rho_basis, const int& nbspline);
    void setup(const UnitCell* ucell,
               const Parallel_Grid& parallel_grid,
               const ModulePW::PW_Basis* rho_basis);

    ModuleBase::ComplexMatrix strucFac;

  private:
    void bspline_sf(int order,
                    const UnitCell* ucell,
                    const Parallel_Grid& parallel_grid,
                    const ModulePW::PW_Basis* rho_basis);
    void bsplinecoef(std::complex<double>* b1,
                     std::complex<double>* b2,
                     std::complex<double>* b3,
                     int nx,
                     int ny,
                     int nz,
                     int order);

    int nbspline = 0;
    const ModulePW::PW_Basis* rho_basis = nullptr;
};

#endif
