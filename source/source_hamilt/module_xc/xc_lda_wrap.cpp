// This file contains wrapper for the LDA functionals
// it includes 3 subroutines:
// 1. xc, which is the wrapper of LDA part
// (i.e. LDA functional and LDA part of GGA functional)
// 2. xc_spin, which is the spin polarized counterpart of xc

#include "xc_functional.h"
#include <stdexcept>

void XC_Functional::xc(
    const double &rho,
    double &exc,
    double &vxc)
{

    double third = 1.0 / 3.0;
    double pi34 = 0.6203504908994e0; // pi34=(3/4pi)^(1/3)
    double rs = 0.0;
    double e = 0.0;
    double v = 0.0;

    exc = 0.00;
    vxc = 0.00;

    rs = pi34 / std::pow(rho, third);

    for (int id : func_id)
    {
        switch (id)
        {
            // Exchange functionals containing slater exchange
            case XC_LDA_X:
            case XC_GGA_X_PBE:
            case XC_GGA_X_PBE_R:
            case XC_GGA_X_PBE_SOL:
            case XC_GGA_X_WC:
            case XC_GGA_X_B88:
            case XC_GGA_X_PW91:
            {
                //  SLA,PBX,rPBX,PBXsol,WC,B88,PW91_X
                XC_Functional::slater(rs, e, v);
                break;
            }

            // Correlation functionals containing PW correlation
            case XC_GGA_C_PBE:
            case XC_GGA_C_PW91:
            case XC_LDA_C_PW:
            case XC_GGA_C_PBE_SOL:
            {
                //   PBC,PW91,PWLDA
                XC_Functional::pw(rs, 0, e, v);
                break;
            }

            // Correlation functionals containing PZ correlation
            case XC_LDA_C_PZ:
            case XC_GGA_C_P86:
            {
                //  PZ,P86
                XC_Functional::pz(rs, 0, e, v);
                break;
            }

            // Correlation functionals containing LYP correlation
            case XC_GGA_C_LYP:
            {
                //  BLYP
                XC_Functional::lyp(rs, e, v);
                break;
            }

            default:
            {
                e = 0.0;
                v = 0.0;
                break;
            }
        }
        exc += e;
        vxc += v;
    }
    return;
}

void XC_Functional::xc_spin(
    const double &rho,
    const double &zeta,
    double &exc,
    double &vxcup,
    double &vxcdw)
{
    static const double small = 1.e-10;
    double e = 0.0;
    double vup = 0.0;
    double vdw = 0.0;
    exc = 0.0;
    vxcup = 0.0;
    vxcdw = 0.0;

    static const double third = 1.0 / 3.0;
    static const double pi34 = 0.62035049089940;
    const double rs = pi34 / pow(rho, third); //wigner_sitz_radius;

    for (int id : func_id)
    {
        switch (id)
        {
            // Exchange functionals containing slater exchange
            case XC_LDA_X:
            case XC_GGA_X_PBE:
            case XC_GGA_X_PBE_R:
            case XC_GGA_X_PBE_SOL:
            case XC_GGA_X_WC:
            case XC_GGA_X_B88:
            case XC_GGA_X_PW91:
            {
                //  SLA,PBX,rPBX,PBXsol,WC,B88,PW91_X
                XC_Functional::slater_spin(rho, zeta, e, vup, vdw);
                break;
            }

            // Correlation functionals containing PZ correlation
            case XC_LDA_C_PZ:
            case XC_GGA_C_P86:
            {
                //  PZ,P86
                XC_Functional::pz_spin(rs, zeta, e, vup, vdw);
                break;
            }

            // Correlation functionals containing PW correlation
            case XC_GGA_C_PBE:
            case XC_GGA_C_PBE_SOL:
            case XC_LDA_C_PW:
            {
                //   PBC,PBCsol
                XC_Functional::pw_spin(rs, zeta, e, vup, vdw);
                break;
            }

            default:
            {
                throw std::domain_error("functional unfinished in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
                break;
            }
        }
        exc += e;
        vxcup += vup;
        vxcdw += vdw;
    }
    return;
}
