#ifndef OPERATOR_LCAO_H
#define OPERATOR_LCAO_H

#include "source_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{

/**
 * Minimal real-space LCAO operator interface used by the H0 executable.
 *
 * The general ABACUS Operator hierarchy owns H(k), wavefunction application,
 * and chained solver state. None of those concepts belongs to the H0-only
 * workflow, so this class deliberately owns only the output H(R) container.
 */
template <typename TK, typename TR>
class OperatorLCAO
{
  public:
    explicit OperatorLCAO(HContainer<TR>* hR_in) : hR(hR_in) {}
    virtual ~OperatorLCAO() {}

    virtual void contributeHR() = 0;

  protected:
    HContainer<TR>* hR = nullptr;
};

} // namespace hamilt

#endif
