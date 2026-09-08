#pragma once
#include <complex>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"

namespace ModuleGint
{

void cal_gint_vl(
    const double* vr_eff,
    hamilt::HContainer<double>* hR);

void cal_gint_vl(
    std::vector<const double*> vr_eff,
    hamilt::HContainer<std::complex<double>>* hR);

} // namespace ModuleGint
