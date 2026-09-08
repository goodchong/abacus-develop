#ifndef READ_PSEUDO_H
#define READ_PSEUDO_H

#include "source_cell/cal_atoms_info.h"
#include "source_cell/unitcell.h"
#include "source_estate/basis_index.h"

namespace elecstate
{

AtomsInfoResult read_pseudo(std::ofstream& ofs,
                            UnitCell& ucell,
                            const std::string& pseudo_dir,
                            const std::string& dft_functional,
                            bool lspinorb,
                            double pseudo_rcut,
                            double soc_lambda,
                            int nspin,
                            int npol,
                            bool two_fermi,
                            double nelec,
                            double nupdown);

void read_cell_pseudopots(const std::string& directory,
                          std::ofstream& log,
                          UnitCell& ucell,
                          const std::string& dft_functional,
                          bool lspinorb,
                          double pseudo_rcut,
                          double soc_lambda);

} // namespace elecstate

#endif
