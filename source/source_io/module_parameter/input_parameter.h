#ifndef INPUT_PARAMETER_H
#define INPUT_PARAMETER_H

#include <string>

// Complete input state for the one-shot LCAO H0 workflow.
struct Input_para
{
    std::string calculation = "get_h0";
    std::string h0_type;
    double h0_sparse_threshold = 1.0e-10;
    int h0_precision = 16;

    std::string suffix = "ABACUS";
    int ntype = 0;
    std::string stru_file = "STRU";
    std::string pseudo_dir;
    std::string orbital_dir;
    std::string read_file_dir = "auto";

    std::string basis_type = "lcao";
    double ecutwfc = 0.0;
    double ecutrho = 0.0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int nbspline = -1;
    double lcao_ecut = 0.0;
    double lcao_dk = 0.01;
    int bx = 0;
    int by = 0;
    int bz = 0;
    int nb2d = 0;

    int nspin = 1;
    bool noncolin = false;
    bool lspinorb = false;
    double soc_lambda = 1.0;
    std::string init_chg = "auto";
    std::string dft_functional = "default";
    double nelec = 0.0;
    double nupdown = 0.0;
    double pseudo_rcut = 15.0;

    double min_dist_coef = 0.2;
};

#endif
