#include "gint_common.h"
#include "source_base/timer.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/tool_quit.h"

#ifdef __MPI
#include "source_base/module_external/blacs_connector.h"
#include <mpi.h>
#endif

namespace ModuleGint
{

template<typename T>
void compose_hr_gint(HContainer<T>& hr_gint)
{
    ModuleBase::TITLE("Gint", "compose_hr_gint");
    ModuleBase::timer::start("Gint", "compose_hr_gint");
    for (int iap = 0; iap < hr_gint.size_atom_pairs(); iap++)
    {
        auto& ap = hr_gint.get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<T>* upper_ap = hr_gint.find_pair(iat2, iat1);
            const hamilt::AtomPair<T>* lower_ap = hr_gint.find_pair(iat1, iat2);
#ifdef __DEBUG
            assert(upper_ap != nullptr);
#endif
            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {
                auto R_index = ap.get_R_index(ir);
                auto upper_mat = upper_ap->find_matrix(-R_index);
                auto lower_mat = lower_ap->find_matrix(R_index);
                for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                    {
                        lower_mat->get_value(icol, irow) = upper_mat->get_value(irow, icol);
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Gint", "compose_hr_gint");
}

template <typename T>
void hr_gint_to_hR(const HContainer<T>& hr_gint, HContainer<T>& hR)
{
    ModuleBase::TITLE("Gint", "hr_gint_to_hR");
    ModuleBase::timer::start("Gint", "hr_gint_to_hR");
#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR.add(hr_gint);
    }
    else
    {
        hamilt::transferSerials2Parallels(hr_gint, &hR);
    }
#else
    hR.add(hr_gint);
#endif
    ModuleBase::timer::end("Gint", "hr_gint_to_hR");
}


void merge_hr_part_to_hR(const std::vector<hamilt::HContainer<double>>& hr_gint_tmp ,
                         hamilt::HContainer<std::complex<double>>* hR,
                         const GintInfo& gint_info){
    ModuleBase::TITLE("Gint_k", "transfer_pvpR");
    ModuleBase::timer::start("Gint_k", "transfer_pvpR");

    const UnitCell* ucell_in = gint_info.get_ucell();
    int mg = hR->get_paraV()->get_global_row_size()/2;
    int ng = hR->get_paraV()->get_global_col_size()/2;
    int nb = hR->get_paraV()->get_block_size()/2;
    hamilt::HContainer<std::complex<double>>* hR_tmp;


#ifdef __MPI
    int blacs_ctxt = hR->get_paraV()->blacs_ctxt;
    std::vector<int> iat2iwt(ucell_in->nat);
    for (int iat = 0; iat < ucell_in->nat; iat++) {
        iat2iwt[iat] = ucell_in->get_iat2iwt()[iat]/2;
    }
    Parallel_Orbitals *pv = new Parallel_Orbitals();
    pv->set(mg, ng, nb, blacs_ctxt);
    pv->set_atomic_trace(iat2iwt.data(), ucell_in->nat, mg);
    auto ijr_info = hR->get_ijr_info();
    hR_tmp = new hamilt::HContainer<std::complex<double>>(pv, nullptr, &ijr_info);
#endif

    //select hr_gint_tmp 
    std::vector<int> first = {0, 1, 1, 0};
    std::vector<int> second= {3, 2, 2, 3};
    //select position in the big matrix
    std::vector<int> row_set = {0, 0, 1, 1};
    std::vector<int> col_set = {0, 1, 0, 1};
    //construct complex matrix
    // Pauli-to-spinor conversion: H = V_0*I + B_x*sigma_x + B_y*sigma_y + B_z*sigma_z
    // sigma_y = [[0,-i],[i,0]], so H_{up,down} = B_x - i*B_y, H_{down,up} = B_x + i*B_y
    // coefficient = clx_i + i*clx_j for each Pauli channel:
    //   is=0 (up,up):   V_0 + B_z   => coeff on B_z = +1  => clx_i=1,  clx_j=0
    //   is=1 (up,down): B_x - i*B_y => coeff on B_y = -i  => clx_i=0,  clx_j=-1
    //   is=2 (down,up): B_x + i*B_y => coeff on B_y = +i  => clx_i=0,  clx_j=+1
    //   is=3 (down,down): -(V_0 - B_z) => coeff on V_0 = -1 => clx_i=-1, clx_j=0
    std::vector<int> clx_i = {1, 0, 0, -1};
    std::vector<int> clx_j = {0, -1, 1, 0};
    for (int is = 0; is < 4; is++){
        if(!PARAM.globalv.domag && (is==1 || is==2)) continue;
#ifdef __MPI
        hR_tmp->set_zero();
#endif
        hamilt::HContainer<std::complex<double>>* hRGint_tmpCd = new hamilt::HContainer<std::complex<double>>(ucell_in->nat);
        hRGint_tmpCd->insert_ijrs( &(gint_info.get_ijr_info()), *(ucell_in));
        hRGint_tmpCd->allocate(nullptr, true);
        hRGint_tmpCd->set_zero();
        for (int iap = 0; iap < hRGint_tmpCd->size_atom_pairs(); iap++)
        {
            auto* ap = &hRGint_tmpCd->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            if (iat1 <= iat2)
            {
                hamilt::AtomPair<std::complex<double>>* upper_ap = ap;
                hamilt::AtomPair<std::complex<double>>* lower_ap = hRGint_tmpCd->find_pair(iat2, iat1);
                const hamilt::AtomPair<double>* ap_nspin1 = hr_gint_tmp [first[is]].find_pair(iat1, iat2);
                const hamilt::AtomPair<double>* ap_nspin2 = hr_gint_tmp [second[is]].find_pair(iat1, iat2);
                for (int ir = 0; ir < upper_ap->get_R_size(); ir++)
                {   
                    const auto R_index = upper_ap->get_R_index(ir);
                    auto upper_mat = upper_ap->find_matrix(R_index);
                    auto mat_nspin1 = ap_nspin1->find_matrix(R_index);
                    auto mat_nspin2 = ap_nspin2->find_matrix(R_index);
                    // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                    for (int irow = 0; irow < mat_nspin1->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < mat_nspin1->get_col_size(); ++icol)
                        {
                            upper_mat->get_value(irow, icol) = mat_nspin1->get_value(irow, icol) 
                            + std::complex<double>(clx_i[is], clx_j[is]) * mat_nspin2->get_value(irow, icol);
                        }
                    }
                    //fill the lower triangle matrix at -R by conjugate transpose of upper at R
                    // This ensures H(-R) = H(R)^dagger, required for Hermiticity of H(k).
                    // For real matrices (is=0,3), conj has no effect.
                    // For complex matrices (is=1,2), conj is essential.
                    if (iat1 < iat2)
                    {
                        auto lower_mat = lower_ap->find_matrix(-R_index);
                        for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                        {
                            for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                            {
                                lower_mat->get_value(icol, irow) = std::conj(upper_mat->get_value(irow, icol));
                            }
                        }
                    }

                }
            }
        }
        // transfer hRGint_tmpCd to parallel hR_tmp
#ifdef __MPI
        hamilt::transferSerials2Parallels( *hRGint_tmpCd, hR_tmp);
#else
        hR_tmp = hRGint_tmpCd;
#endif
        // merge hR_tmp to hR
        for (int iap = 0; iap < hR->size_atom_pairs(); iap++)
        {
            auto* ap = &hR->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            auto* ap_nspin = hR_tmp ->find_pair(iat1, iat2);
            if (ap_nspin == nullptr)
            {
                continue;
            }
            for (int ir = 0; ir < ap->get_R_size(); ir++)
            {   
                const auto R_index = ap->get_R_index(ir);
                auto upper_mat = ap->find_matrix(R_index);
                auto mat_nspin = ap_nspin->find_matrix(R_index);
                if (mat_nspin == nullptr)
                {
                    continue;
                }
                // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                for (int irow = 0; irow < mat_nspin->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < mat_nspin->get_col_size(); ++icol)
                    {
                        upper_mat->get_value(2*irow+row_set[is], 2*icol+col_set[is]) = 
                        mat_nspin->get_value(irow, icol);
                    }
                }
            }
        }
        delete hRGint_tmpCd;
    }
#ifdef __MPI
    delete hR_tmp;
    delete pv;
#endif
    ModuleBase::timer::end("Gint_k", "transfer_pvpR");
    return;
}




template void compose_hr_gint(HContainer<double>& hr_gint);
template void compose_hr_gint(HContainer<float>& hr_gint);
template void hr_gint_to_hR(
    const HContainer<double>& hr_gint,
    HContainer<double>& hR);
template void hr_gint_to_hR(
    const HContainer<std::complex<double>>& hr_gint,
    HContainer<std::complex<double>>& hR);
}
