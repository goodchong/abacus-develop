#pragma once
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_gint/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    template<typename T>
    void compose_hr_gint(HContainer<T>& hr_gint);

    template <typename T>
    void hr_gint_to_hR(const HContainer<T>& hr_gint, HContainer<T>& hR);
    // for nspin=4 case
    void merge_hr_part_to_hR(const std::vector<hamilt::HContainer<double>>& hr_gint_tmp ,
                         hamilt::HContainer<std::complex<double>>* hR,
                         const GintInfo& gint_info);

}
