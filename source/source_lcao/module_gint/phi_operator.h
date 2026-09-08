#pragma once

#include "big_grid.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace ModuleGint
{

/** Operations required to project a local potential onto an H(R) container. */
class PhiOperator
{
  public:
    enum class TriPart
    {
        Upper,
        Lower,
        Full
    };

    void set_bgrid(std::shared_ptr<const BigGrid> biggrid);

    int get_rows() const { return rows_; }
    int get_cols() const { return cols_; }

    template <typename T>
    void set_phi(T* phi) const;

    template <typename T>
    void phi_mul_vldr3(const T* vl, const T dr3, const T* phi, T* result) const;

    template <typename Tin>
    void phi_mul_phi(const Tin* phi_i,
                     const Tin* phi_j,
                     HContainer<double>& hr,
                     const TriPart part) const;

  private:
    void init_atom_pair_idx_();

    const std::pair<int, int>& get_atom_pair_start_end_idx_(int a, int b) const
    {
        const int x = std::min(a, b);
        const int y = std::abs(a - b);
        return atom_pair_range_[(2 * biggrid_->get_atoms_num() - x + 1) * x / 2 + y];
    }

    bool is_atom_on_mgrid(int atom_idx, int mgrid_idx) const
    {
        return is_atom_on_mgrid_[atom_idx * rows_ + mgrid_idx];
    }

    int rows_ = 0;
    int cols_ = 0;
    std::vector<int> mgrid_lidx_;
    std::shared_ptr<const BigGrid> biggrid_;
    std::vector<std::vector<Vec3d>> atom_rcoords_;
    std::vector<bool> is_atom_on_mgrid_;
    std::vector<int> atoms_startidx_;
    std::vector<int> atoms_phi_len_;
    std::vector<std::pair<int, int>> atom_pair_range_;
};

} // namespace ModuleGint

#include "phi_operator.hpp"
