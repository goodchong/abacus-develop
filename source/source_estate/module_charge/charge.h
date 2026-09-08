#ifndef CHARGE_H
#define CHARGE_H

#include "source_base/complexmatrix.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_basis/module_pw/pw_basis.h"
// #include "source_estate/fp_energy.h"
#include "source_pw/module_pwdft/parallel_grid.h"

//a forward declaration of UnitCell
class UnitCell;

// Electron Charge Density
class Charge
{

  public:

    Charge();
    ~Charge();

    //==========================================================
    // MEMBER VARIABLES :
    // init_chg : "atomic" or "file"
    // NAME : total number of electrons
    // NAME : rho (nspin,ncxyz), the charge density in real space
    // NAME : rhog, charge density in G space
    // NAME : rho_core [nrxx], the core charge in real space
    // NAME : rhog_core [ngm], the core charge in reciprocal space
    //==========================================================

    double **rho = nullptr;
    std::complex<double> **rhog = nullptr;
    const Parallel_Grid* pgrid = nullptr;

  private:

    //temporary
    double *_space_rho = nullptr; 
    std::complex<double> *_space_rhog = nullptr;

  public:

    double *rho_core = nullptr;
    std::complex<double> *rhog_core = nullptr;

    void set_rhopw(ModulePW::PW_Basis* rhopw_in);

    /**
     * @brief Init charge density from file or atomic pseudo-wave-functions
     *
     * @param eferm_iout [out] fermi energy to be initialized
     * @param ucell [in] unit cell
     * @param strucFac [in] structure factor
     */
    void init_rho(const UnitCell& ucell,
                  const Parallel_Grid& pgrid,
                  const ModuleBase::ComplexMatrix& strucFac);

    void allocate(int nspin_in);

    void atomic_rho(const int spin_number_need,
                    const double& omega,
                    double** rho_in,
                    const ModuleBase::ComplexMatrix& strucFac,
                    const UnitCell& ucell) const;

    void set_rho_core(const UnitCell& ucell,
                      const ModuleBase::ComplexMatrix& structure_factor, 
                      const bool* numeric);

    void renormalize_rho();

    double sum_rho() const;

	// for non-linear core correction
    void non_linear_core_correction
    (
        const bool &numeric,
        const double omega,
        const double tpiba2,
        const int mesh,
        const double *r,
        const double *rab,
        const double *rhoc,
        double *rhocg
    ) const;

	double cal_rho2ne(const double *rho_in) const;

    void check_rho(); // to check whether the charge density is normal

    void set_omega(double* omega_in){this->omega_ = omega_in;};

    // mohan add 2021-02-20
    int nrxx=0; // number of r vectors in this processor
    int ngmc=0; // number of g vectors in this processor
    int nspin=0; // number of spins
    ModulePW::PW_Basis* rhopw = nullptr;// When double_grid is used, rhopw = rhodpw (dense grid)

  private:

    void destroy();    // free arrays  liuyu 2023-03-12

    double* omega_ = nullptr; // omega for non-linear core correction

    bool allocate_rho = false;
    
};

#endif // charge
