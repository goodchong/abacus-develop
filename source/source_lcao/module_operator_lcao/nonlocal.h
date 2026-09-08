#ifndef NONLOCALNEW_H
#define NONLOCALNEW_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <unordered_map>
#include <vector>

namespace hamilt
{

#ifndef __NONLOCALNEWTEMPLATE
#define __NONLOCALNEWTEMPLATE

template <class T>
class Nonlocal : public T
{
};

#endif

/// Nonlocal class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the non-local pseudopotential matrix in real space.
/// HR = <psi_{mu, 0}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, R}>
template <typename TK, typename TR>
class Nonlocal<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    Nonlocal<OperatorLCAO<TK, TR>>(hamilt::HContainer<TR>* hR_in,
                                      const UnitCell* ucell_in,
                                      const std::vector<double>& orb_cutoff,
                                      const Grid_Driver* GridD_in,
                                      const TwoCenterIntegrator* intor);
    ~Nonlocal<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|beta_p1>D_{p1, p2}<beta_p2|phi_{\nu, R}>
     */
    void contributeHR() override;

  private:
    const UnitCell* ucell = nullptr;

    std::vector<double> orb_cutoff_;

    hamilt::HContainer<TR>* HR_fixed = nullptr;

    // the following variable is introduced temporarily during LCAO refactoring
    const TwoCenterIntegrator* intor_ = nullptr;

    bool allocated = false;

    bool HR_fixed_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const Grid_Driver* GridD_in);

    /**
     * @brief calculate the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in HR and calculate the non-local pseudopotential matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const int& T0,
                    const Parallel_Orbitals* paraV,
                    const std::unordered_map<int, std::vector<double>>& nlm1_all,
                    const std::unordered_map<int, std::vector<double>>& nlm2_all,
                    TR* data_pointer);

    const Grid_Driver* gridD = nullptr;

    std::vector<AdjacentAtomInfo> adjs_all;
};

} // namespace hamilt
#endif
