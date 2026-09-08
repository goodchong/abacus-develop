#ifndef EKINETIC_H
#define EKINETIC_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include <vector>

namespace hamilt
{

#ifndef __EKINETICTEMPLATE
#define __EKINETICTEMPLATE

template <class T>
class EKinetic : public T
{
};

#endif

/// EKinetic class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the electronic kinetic matrix in real space.
/// HR = <psi_{mu, 0}|-\Nabla^2|psi_{nu, R}>
template <typename TK, typename TR>
class EKinetic<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    /**
     * @brief Construct a new EKinetic object
     */
    EKinetic<OperatorLCAO<TK, TR>>(HContainer<TR>* hR_in,
                                      const UnitCell* ucell_in,
                                      const std::vector<double>& orb_cutoff,
                                      const Grid_Driver* GridD_in,
                                      const TwoCenterIntegrator* intor);

    /**
     * @brief Destroy the EKinetic object
     */
    ~EKinetic<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|-\Nabla^2|phi_{\nu, R}>
     */
    void contributeHR() override;

  private:
    const UnitCell* ucell = nullptr;
    std::vector<double> orb_cutoff_;

    hamilt::HContainer<TR>* HR_fixed = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    const Grid_Driver* gridD = nullptr;

    bool allocated = false;

    bool HR_fixed_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const Grid_Driver* GridD_in);

    /**
     * @brief calculate the electronic kinetic matrix with specific <I,J,R> atom-pairs
     * use the adjs_all to calculate the HR matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const ModuleBase::Vector3<double>& dtau,
                    TR* data_pointer);

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_all;
};

} // namespace hamilt
#endif
