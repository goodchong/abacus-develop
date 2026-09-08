#ifndef TWO_CENTER_BUNDLE_H
#define TWO_CENTER_BUNDLE_H

#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_nao/two_center_integrator.h"

#include <memory>
#include <string>

class TwoCenterBundle
{
  public:
    TwoCenterBundle() = default;
    ~TwoCenterBundle() = default;
    TwoCenterBundle& operator=(TwoCenterBundle&&) = default;

    // NOTE: some variables might be set only on RANK-0
    void build_orb(int ntype, const std::string* file_orb0, const std::string& orbital_dir);
    void build_beta(int ntype, Numerical_Nonlocal* nl);

    void tabulate();

    /**
     * @brief Overwrites the content of a LCAO_Orbitals object (e.g. ORB)
     * with the current object.
     *
     * This function provides an interface to the corresponding object in the old module_ao.
     */
    void to_LCAO_Orbitals(LCAO_Orbitals& orb,
                          const double lcao_ecut,
                          const double lcao_dk,
                          const bool out_element_info) const;

    std::unique_ptr<TwoCenterIntegrator> kinetic_orb;
    std::unique_ptr<TwoCenterIntegrator> overlap_orb_beta;

    std::unique_ptr<RadialCollection> orb_;
    std::unique_ptr<RadialCollection> beta_;
};

#endif
