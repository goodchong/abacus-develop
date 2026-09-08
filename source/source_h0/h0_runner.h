#ifndef H0_RUNNER_H
#define H0_RUNNER_H

#include "source_cell/unitcell.h"

#include <memory>

struct Input_para;

namespace ModuleH0
{

/**
 * @brief One-shot LCAO initial-Hamiltonian builder.
 *
 * The runner owns only the state needed to form H(R), writes the requested
 * matrix, and returns without any electronic solver or iterative state.
 */
class H0Runner final
{
  public:
    H0Runner();
    ~H0Runner();

    void initialize(UnitCell& ucell, const Input_para& inp);
    void run(UnitCell& ucell);
    void finalize();

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ModuleH0

#endif
