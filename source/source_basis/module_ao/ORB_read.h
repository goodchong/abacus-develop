#ifndef LCAO_ORBITALS_H
#define LCAO_ORBITALS_H

#include "ORB_atomic.h"

#include <vector>

class LCAO_Orbitals
{
  public:
    LCAO_Orbitals();
    ~LCAO_Orbitals();

    LCAO_Orbitals(const LCAO_Orbitals&) = delete;
    LCAO_Orbitals& operator=(const LCAO_Orbitals&) = delete;

    const int& get_kmesh() const { return kmesh; }
    const double& get_dk() const { return dk; }
    const int& get_ntype() const { return ntype; }
    const double& get_dr_uniform() const { return dr_uniform; }
    const double& get_rcutmax_Phi() const { return rcutmax_Phi; }

    std::vector<double> cutoffs() const;

    Numerical_Orbital* Phi = nullptr;

  private:
    int ntype = 0;
    int kmesh = 0;
    double dk = 0.01;
    double dr_uniform = 0.001;
    double rcutmax_Phi = 0.0;

    friend class TwoCenterBundle;
};

#endif
