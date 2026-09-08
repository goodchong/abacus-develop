#include "phi_operator.h"

namespace ModuleGint
{

void PhiOperator::set_bgrid(std::shared_ptr<const BigGrid> biggrid)
{
    biggrid_ = biggrid;
    rows_ = biggrid_->get_mgrids_num();
    cols_ = biggrid_->get_phi_len();

    biggrid_->set_atoms_startidx(atoms_startidx_);
    biggrid_->set_atoms_phi_len(atoms_phi_len_);
    biggrid_->set_mgrids_local_idx(mgrid_lidx_);

    const int atoms_num = biggrid_->get_atoms_num();
    atom_rcoords_.resize(atoms_num);
    is_atom_on_mgrid_.resize(rows_ * atoms_num);
    for (int i = 0; i < atoms_num; ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        biggrid_->set_atom_relative_coords(atom, atom_rcoords_[i]);
        for (int j = 0; j < rows_; ++j)
        {
            is_atom_on_mgrid_[i * rows_ + j] = atom_rcoords_[i][j].norm() <= atom->get_rcut();
        }
    }
    init_atom_pair_idx_();
}

void PhiOperator::init_atom_pair_idx_()
{
    const int atoms_num = biggrid_->get_atoms_num();
    const int mgrids_num = biggrid_->get_mgrids_num();
    atom_pair_range_.resize(atoms_num * (atoms_num + 1) / 2);
    int atom_pair_idx = 0;
    for (int i = 0; i < atoms_num; ++i)
    {
        for (int j = i; j < atoms_num; ++j)
        {
            int start_idx = mgrids_num;
            int end_idx = -1;
            for (int mgrid_idx = 0; mgrid_idx < mgrids_num; ++mgrid_idx)
            {
                if (is_atom_on_mgrid(i, mgrid_idx) && is_atom_on_mgrid(j, mgrid_idx))
                {
                    start_idx = mgrid_idx;
                    break;
                }
            }
            for (int mgrid_idx = mgrids_num - 1; mgrid_idx >= 0; --mgrid_idx)
            {
                if (is_atom_on_mgrid(i, mgrid_idx) && is_atom_on_mgrid(j, mgrid_idx))
                {
                    end_idx = mgrid_idx;
                    break;
                }
            }
            atom_pair_range_[atom_pair_idx++] = std::make_pair(start_idx, end_idx);
        }
    }
}

} // namespace ModuleGint
