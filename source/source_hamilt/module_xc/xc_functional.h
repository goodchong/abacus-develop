//==========================================================
// AUTHOR : mohan
// DATE : 2009-04-04
//==========================================================
#ifndef XC_FUNCTIONAL_H
#define XC_FUNCTIONAL_H

#include "xc_ids.h"
#include "source_base/macros.h"
#include "source_base/global_function.h"
#include "source_base/vector3.h"
#include "source_base/matrix.h"
#include "source_estate/module_charge/charge.h"
#include "source_cell/unitcell.h"

class XC_Functional
{
    public:

    XC_Functional();
    ~XC_Functional();

//-------------------
// subroutines, grouped according to the file they are in:
//-------------------

//-------------------
//  xc_pot.cpp
//-------------------

    // Compute the built-in LDA/GGA exchange-correlation potential.
    static std::tuple<double, double, ModuleBase::matrix> v_xc(
        const int &nrxx, // number of real-space grid
        const Charge* const chr,
        const UnitCell *ucell, // charge density
        const int nspin,
        const bool domag,
        const bool domag_z);

//-------------------
//  xc_functional.cpp
//-------------------

    static int get_func_type()
    {
        return func_type;
    };

    static void set_xc_type(const std::string xc_func_in);

    static bool get_ked_flag()
    {
        return false;
    };

    static std::string output_info();

    private:

    static std::vector<int> func_id;
    static int func_type; // 1 = LDA, 2 = GGA

//-------------------
//  xc_lda_wrap.cpp
//-------------------

// This file contains wrappers for the built-in LDA functionals:
// 1. xc, which is the wrapper of LDA part
// (i.e. LDA functional and LDA part of GGA functional)
// 2. xc_spin, which is the spin polarized counterpart of xc

    public:

    // LDA
    static void xc(const double &rho, double &exc, double &vxc);

    // LSDA
    static void xc_spin(
        const double &rho,
        const double &zeta,
        double &exc,
        double &vxcup,
        double &vxcdw);

//-------------------
//  xc_gga_wrap.cpp
//-------------------

// This file contains wrapper for the GGA functionals
// it includes 4 subroutines:
// 1. gcxc, which is the wrapper for gradient correction part
// 2. gcx_spin, spin polarized, exchange only
// 3. gcc_spin, spin polarized, correlation only

    // GGA
    static void gcxc(
        const double &rho,
        const double &grho,
        double &sxc,
        double &v1xc,
        double &v2xc);

    // spin polarized GGA
    static void gcx_spin(
        double rhoup,
        double rhodw,
        double grhoup2,
        double grhodw2,
        double &sx,
        double &v1xup,
        double &v1xdw,
        double &v2xup,
        double &v2xdw);

    static void gcc_spin(
        double rho,
        double &zeta,
        double grho,
        double &sc,
        double &v1cup,
        double &v1cdw,
        double &v2c);

//-------------------
//  xc_grad.cpp
//-------------------

// This file contains the gradient calculations used by built-in GGA:
// 1. gradcorr, which calculates gradient correction
// 2. grad_rho, which calculates gradient of density
// 3. grad_dot, which calculates divergence of something
// 4. noncolin_rho, which diagonalizes the spin density matrix
//  and gives the spin up and spin down components of the charge.

    static void gradcorr(
        double& etxc,
        double& vtxc,
        ModuleBase::matrix& v,
        const Charge* const chr,
        ModulePW::PW_Basis* rhopw,
        const UnitCell* ucell,
        const int nspin,
        const bool domag,
        const bool domag_z);

    static void grad_rho(
        const std::complex<double>* rhog,
        ModuleBase::Vector3<double>* gdr,
        const ModulePW::PW_Basis* rho_basis,
        const double tpiba);

    static void grad_dot(
        const ModuleBase::Vector3<double>* h,
        double* dh,
        const ModulePW::PW_Basis* rho_basis,
        const double tpiba);

    static void noncolin_rho(
        double* rhoout1,
        double* rhoout2,
        double* seg,
        const double* const* const rho,
        const int nrxx,
        const double* ux_,
        const bool lsign_);

    //-------------------
    //  xc_lda_exch.cpp
    //-------------------

    // This file contains realization of LDA exchange functionals
    // Spin unpolarized ones:
    //  1. slater: ordinary Slater exchange with alpha=2/3
    //  2. slater1: Slater exchange with alpha=1
    //  3. slater_rxc : Slater exchange with alpha=2/3 and Relativistic exchange
    // And their spin polarized counterparts:
    //  1. slater_spin
    //  2. slater1_spin
    //  3. slater_rxc_spin

    // For LDA exchange energy
    static void slater(const double &rs, double &ex, double &vx);
    static void slater1(const double &rs, double &ex, double &vx);
    static void slater_rxc(const double &rs, double &ex, double &vx);

    // For LSDA exchange energy
    static void slater_spin(
        const double &rho,
        const double &zeta,
        double &ex,
        double &vxup,
        double &vxdw);
    static void slater1_spin(
        const double &rho,
        const double &zeta,
        double &ex,
        double &vxup,
        double &vxdw);
    static void slater_rxc_spin(
        const double &rho,
        const double &z,
        double &ex,
        double &vxup,
        double &vxdw);

//-------------------
//  xc_lda_corr.cpp
//-------------------

// This file contains realization of LDA correlation functionals
// Spin unpolarized ones:
//  1. pw : Perdew-Wang LDA correlation
//  2. pz : Perdew-Zunger LDA correlation
//  3. lyp : Lee-Yang-Parr correlation
//  4. vwn : Vosko-Wilk-Nusair LDA correlation
//  5. wigner : Wigner
//  6. hl : Hedin-Lunqvist
//  7. gl : Gunnarson-Lunqvist
// And some of their spin polarized counterparts:
//  1. pw_spin
//  2. pz_spin, which calls pz_polarized

    // For LDA correlation energy
    static void pw(const double &rs, const int &iflag, double &ec, double &vc);
    static void pz(const double &rs, const int &iflag, double &ec, double &vc);
    static void lyp(const double &rs, double &ec, double &vc);
    static void vwn(const double &rs, double &ec, double &vc);
    static void wigner(const double &rs, double &ec, double &vc);
    static void hl(const double &rs, double &ec, double &vc);
    static void gl(const double &rs, double &ec, double &vc);

    // For LSDA correlation energy
    static void pw_spin(
        const double &rs,
        const double &zeta,
        double &ec,
        double &vcup,
        double &vcdw);

    static void pz_spin(
        const double &rs,
        const double &zeta,
        double &ec,
        double &vcup,
        double &vcdw);

    static void pz_polarized(const double &rs, double &ec, double &vc);

//-------------------
//  xc_gga_exch.cpp
//-------------------

// This file contains realizations of gradient correction to exchange part
// Spin unpolarized ones:
//  1. becke88 : Becke88 exchange
//  2. ggax : PW91 exchange
//  3. pbex : PBE exchange (and revPBE)
//  4. optx : OPTX, Handy et al.
//  5. wcx : Wu-Cohen exchange
// And some of their spin polarized counterparts:
//  1. becke88_spin

    static void becke88(
        const double &rho,
        const double &grho,
        double &sx,
        double &v1x,
        double &v2x);

    static void ggax(
        const double &rho,
        const double &grho,
        double &sx,
        double &v1x,
        double &v2x);

    static void pbex(
        const double &rho,
        const double &grho,
        const int &iflag,
        double &sx,
        double &v1x,
        double &v2x);

    static void optx(const double rho, const double grho, double &sx, double &v1x, double &v2x);

    static void wcx(
        const double &rho,
        const double &grho,
        double &sx,
        double &v1x,
        double &v2x);

    static void becke88_spin(
        double rho,
        double grho,
        double &sx,
        double &v1x,
        double &v2x);

//-------------------
//  xc_gga_corr.cpp
//-------------------

// This file contains realizations of gradient correction to correlation part
// Spin unpolarized ones:
//  1. perdew86 : P86
//  2. ggac : PW91
//  3. pbec
//  4. glyp
// And some of their spin polarized counterparts:
//  1. perdew86_spin
//  2. ggac_spin
//  3. pbec_spin

    static void perdew86(const double rho, const double grho, double &sc, double &v1c, double &v2c);

    static void ggac(
        const double &rho,
        const double &grho,
        double &sc,
        double &v1c,
        double &v2c);

    static void pbec(
        const double &rho,
        const double &grho,
        const int &flag,
        double &sc,
        double &v1c,
        double &v2c);

    static void glyp(
        const double &rho,
        const double &grho,
        double &sc,
        double &v1c,
        double &v2c);

    static void perdew86_spin(
        double rho,
        double zeta,
        double grho,
        double &sc,
        double &v1cup,
        double &v1cdw,
        double &v2c);

    //static void ggac_spin(double rho, double zeta, double grho, double &sc,
    //  double &v1cup, double &v1cdw, double &v2c);

    static void pbec_spin(
        double rho,
        double zeta,
        double grho,
        const int &flag,
        double &sc,
        double &v1cup,
        double &v1cdw,
        double &v2c);

//-------------------
//  xc_funct_hcth.cpp
//-------------------
// This file contains realizations of the HCTH GGA functional
// hcth calls pwcorr

    static void hcth(const double rho, const double grho, double &sx, double &v1x, double &v2x);

    static void pwcorr(const double r, const double c[], double &g, double &dg);

};

#endif //XC_FUNCTION_H
