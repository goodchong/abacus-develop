#pragma once

#include "source_base/global_function.h"
#include "source_base/module_external/blas_connector.h"

namespace ModuleGint
{

template <typename T>
void PhiOperator::set_phi(T* phi) const
{
    for (int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom = biggrid_->get_atom(i);
        atom->set_phi(atom_rcoords_[i], cols_, phi);
        phi += atom->get_nw();
    }
}

template <typename T>
void PhiOperator::phi_mul_vldr3(const T* vl, const T dr3, const T* phi, T* result) const
{
    int idx = 0;
    for (int i = 0; i < rows_; ++i)
    {
        const T vldr3_mgrid = vl[mgrid_lidx_[i]] * dr3;
        for (int j = 0; j < cols_; ++j)
        {
            result[idx] = phi[idx] * vldr3_mgrid;
            ++idx;
        }
    }
}

template <typename Tin>
void PhiOperator::phi_mul_phi(const Tin* phi_i,
                              const Tin* phi_j,
                              HContainer<double>& hr,
                              const TriPart part) const
{
    std::vector<Tin> tmp_hr;
    for (int i = 0; i < biggrid_->get_atoms_num(); ++i)
    {
        const auto atom_i = biggrid_->get_atom(i);
        const auto& r_i = atom_i->get_R();
        const int iat_i = atom_i->get_iat();
        const int n_i = atoms_phi_len_[i];

        for (int j = 0; j < biggrid_->get_atoms_num(); ++j)
        {
            const auto atom_j = biggrid_->get_atom(j);
            const auto& r_j = atom_j->get_R();
            const int iat_j = atom_j->get_iat();
            const int n_j = atoms_phi_len_[j];

            if ((part == TriPart::Upper && iat_i > iat_j)
                || (part == TriPart::Lower && iat_i < iat_j))
            {
                continue;
            }

            auto* result = hr.find_matrix(iat_i, iat_j, r_i - r_j);
            if (result == nullptr)
            {
                continue;
            }

            const int start_idx = get_atom_pair_start_end_idx_(i, j).first;
            const int len = get_atom_pair_start_end_idx_(i, j).second - start_idx + 1;
            if (len <= 0)
            {
                continue;
            }

            tmp_hr.assign(n_i * n_j, Tin(0));
            const Tin alpha = 1;
            const Tin beta = 1;
            BlasConnector::gemm('T',
                                'N',
                                n_i,
                                n_j,
                                len,
                                alpha,
                                phi_i + start_idx * cols_ + atoms_startidx_[i],
                                cols_,
                                phi_j + start_idx * cols_ + atoms_startidx_[j],
                                cols_,
                                beta,
                                tmp_hr.data(),
                                n_j,
                                base_device::AbacusDevice_t::CpuDevice);
            result->add_array_ts(tmp_hr.data());
        }
    }
}

} // namespace ModuleGint
