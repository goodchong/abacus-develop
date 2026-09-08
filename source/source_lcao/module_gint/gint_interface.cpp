#include "gint_interface.h"

#include "gint_vl.h"
#include "gint_vl_nspin4.h"

namespace ModuleGint
{

void cal_gint_vl(const double* vr_eff, hamilt::HContainer<double>* hR)
{
    Gint_vl gint_vl(vr_eff, hR);
    gint_vl.cal_gint();
}

void cal_gint_vl(std::vector<const double*> vr_eff, hamilt::HContainer<std::complex<double>>* hR)
{
    Gint_vl_nspin4 gint_vl_nspin4(vr_eff, hR);
    gint_vl_nspin4.cal_gint();
}

} // namespace ModuleGint
