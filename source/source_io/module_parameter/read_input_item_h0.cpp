#include "read_input.h"
#include "read_input_tool.h"

#include "source_base/global_function.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

namespace ModuleIO
{

namespace
{

std::string lowercase(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

void require_positive(const char* label, const double value)
{
    if (value <= 0.0)
    {
        ModuleBase::WARNING_QUIT("ReadInput", std::string(label) + " must be positive.");
    }
}

} // namespace

void ReadInput::item_h0()
{
    {
        Input_Item item("calculation");
        item.annotation = "one-shot H0 matrix construction";
        item.category = "H0 workflow";
        item.type = "String";
        item.description = "The H0-only executable accepts get_h0 and exits after writing H(R).";
        item.default_value = "get_h0";
        read_sync_string(input.calculation);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.calculation != "get_h0")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "calculation must be get_h0.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("h0_type");
        item.annotation = "core or full";
        item.category = "H0 workflow";
        item.type = "String";
        item.description = "Required. core constructs T + V_nl + V_loc. full additionally constructs V_H[rho0] + V_xc[rho0] from the initial density, without SCF.";
        item.default_value = "Required";
        read_sync_string(input.h0_type);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read() || (para.input.h0_type != "core" && para.input.h0_type != "full"))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "h0_type must be explicitly set to core or full.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("h0_sparse_threshold");
        item.annotation = "absolute CSR filtering threshold";
        item.category = "H0 workflow";
        item.type = "Real";
        item.description = "Write a matrix entry only when abs(value) is strictly greater than this threshold.";
        item.default_value = "1e-10";
        read_sync_double(input.h0_sparse_threshold);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.h0_sparse_threshold < 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "h0_sparse_threshold must be non-negative.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("h0_precision");
        item.annotation = "significant digits in text CSR values";
        item.category = "H0 workflow";
        item.type = "Integer";
        item.description = "Number of significant digits written for matrix values; valid range is 1 through 17.";
        item.default_value = "16";
        read_sync_int(input.h0_precision);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.h0_precision < 1 || para.input.h0_precision > 17)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "h0_precision must be between 1 and 17.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("suffix");
        item.annotation = "output directory suffix";
        item.category = "Files";
        item.type = "String";
        item.description = "Write results beneath OUT.<suffix>.";
        item.default_value = "ABACUS";
        read_sync_string(input.suffix);
        this->add_item(item);
    }
    {
        Input_Item item("ntype");
        item.annotation = "number of atomic species";
        item.category = "Files";
        item.type = "Integer";
        item.description = "Number of atomic species. Zero detects the value from STRU.";
        item.default_value = "0";
        read_sync_int(input.ntype);
        this->add_item(item);
    }
    {
        Input_Item item("stru_file");
        item.annotation = "structure filename";
        item.category = "Files";
        item.type = "String";
        item.description = "Structure file containing species, pseudopotential and orbital filenames, lattice, and atomic positions.";
        item.default_value = "STRU";
        read_sync_string(input.stru_file);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_dir");
        item.annotation = "pseudopotential directory";
        item.category = "Files";
        item.type = "String";
        item.description = "Directory prepended to pseudopotential filenames in STRU.";
        item.default_value = "./";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.pseudo_dir = item.get_size() == 0 ? "" : to_dir(strvalue);
        };
        sync_string(input.pseudo_dir);
        this->add_item(item);
    }
    {
        Input_Item item("orbital_dir");
        item.annotation = "numerical-orbital directory";
        item.category = "Files";
        item.type = "String";
        item.description = "Directory prepended to numerical atomic orbital filenames in STRU.";
        item.default_value = "./";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.orbital_dir = item.get_size() == 0 ? "" : to_dir(strvalue);
        };
        sync_string(input.orbital_dir);
        this->add_item(item);
    }
    {
        Input_Item item("read_file_dir");
        item.annotation = "initial-density input directory";
        item.category = "Files";
        item.type = "String";
        item.description = "Directory searched for ABACUS restart density and spin-resolved cube files when init_chg is file or auto.";
        item.default_value = "OUT.<suffix>";
        read_sync_string(input.read_file_dir);
        item.reset_value = [](const Input_Item&, Parameter& para) {
            if (para.input.read_file_dir == "auto")
            {
                para.input.read_file_dir = "OUT." + para.input.suffix;
            }
            para.input.read_file_dir = to_dir(para.input.read_file_dir);
        };
        this->add_item(item);
    }
    {
        Input_Item item("basis_type");
        item.annotation = "localized atomic orbitals";
        item.category = "Basis";
        item.type = "String";
        item.description = "Only lcao is supported.";
        item.default_value = "lcao";
        read_sync_string(input.basis_type);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.basis_type != "lcao")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "basis_type must be lcao.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("ecutwfc");
        item.annotation = "FFT and radial-table cutoff";
        item.category = "Basis";
        item.type = "Real";
        item.unit = "Ry";
        item.description = "Plane-wave cutoff used for real-space integration grids and LCAO radial tables.";
        item.default_value = "100";
        read_sync_double(input.ecutwfc);
        item.reset_value = [](const Input_Item&, Parameter& para) {
            if (para.input.ecutwfc == 0.0)
            {
                para.input.ecutwfc = para.input.ecutrho > 0.0 ? para.input.ecutrho / 4.0 : 100.0;
            }
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            require_positive("ecutwfc", para.input.ecutwfc);
        };
        this->add_item(item);
    }
    {
        Input_Item item("ecutrho");
        item.annotation = "charge and potential cutoff";
        item.category = "Basis";
        item.type = "Real";
        item.unit = "Ry";
        item.description = "Charge/potential grid cutoff. The LCAO H0 workflow requires ecutrho = 4*ecutwfc unless nx, ny, and nz are supplied explicitly.";
        item.default_value = "4*ecutwfc";
        read_sync_double(input.ecutrho);
        item.reset_value = [](const Input_Item&, Parameter& para) {
            if (para.input.ecutrho <= 0.0)
            {
                para.input.ecutrho = 4.0 * para.input.ecutwfc;
            }
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const bool automatic_grid = para.input.nx * para.input.ny * para.input.nz == 0;
            if (automatic_grid
                && std::abs(para.input.ecutrho / para.input.ecutwfc - 4.0) > 1.0e-8)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ecutrho must equal 4*ecutwfc for the LCAO H0 workflow.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nx");
        item.annotation = "FFT grid size along x";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Explicit FFT grid dimension; nx, ny, and nz must be supplied together.";
        item.default_value = "0";
        read_sync_int(input.nx);
        this->add_item(item);
    }
    {
        Input_Item item("ny");
        item.annotation = "FFT grid size along y";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Explicit FFT grid dimension; nx, ny, and nz must be supplied together.";
        item.default_value = "0";
        read_sync_int(input.ny);
        this->add_item(item);
    }
    {
        Input_Item item("nz");
        item.annotation = "FFT grid size along z";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Explicit FFT grid dimension; nx, ny, and nz must be supplied together.";
        item.default_value = "0";
        read_sync_int(input.nz);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const bool none = para.input.nx == 0 && para.input.ny == 0 && para.input.nz == 0;
            const bool all = para.input.nx > 0 && para.input.ny > 0 && para.input.nz > 0;
            if (!none && !all)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, and nz must be all zero or all positive.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nbspline");
        item.annotation = "structure-factor B-spline order";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Non-negative values enable cardinal B-spline structure factors; -1 disables them.";
        item.default_value = "-1";
        read_sync_int(input.nbspline);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.nbspline < -1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nbspline must be -1 or non-negative.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("lcao_ecut");
        item.annotation = "two-center table cutoff";
        item.category = "Basis";
        item.type = "Real";
        item.unit = "Ry";
        item.description = "Upper reciprocal-space cutoff for LCAO two-center integral tables.";
        item.default_value = "ecutwfc";
        read_sync_double(input.lcao_ecut);
        item.reset_value = [](const Input_Item&, Parameter& para) {
            if (para.input.lcao_ecut == 0.0)
            {
                para.input.lcao_ecut = para.input.ecutwfc;
            }
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            require_positive("lcao_ecut", para.input.lcao_ecut);
        };
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dk");
        item.annotation = "two-center reciprocal spacing";
        item.category = "Basis";
        item.type = "Real";
        item.description = "Reciprocal-space spacing for two-center integral tables.";
        item.default_value = "0.01";
        read_sync_double(input.lcao_dk);
        item.check_value = [](const Input_Item&, const Parameter& para) { require_positive("lcao_dk", para.input.lcao_dk); };
        this->add_item(item);
    }
    {
        Input_Item item("bx");
        item.annotation = "Gint block size along x";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Real-space integration block size; zero selects it automatically.";
        item.default_value = "0";
        read_sync_int(input.bx);
        this->add_item(item);
    }
    {
        Input_Item item("by");
        item.annotation = "Gint block size along y";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Real-space integration block size; zero selects it automatically.";
        item.default_value = "0";
        read_sync_int(input.by);
        this->add_item(item);
    }
    {
        Input_Item item("bz");
        item.annotation = "Gint block size along z";
        item.category = "Basis";
        item.type = "Integer";
        item.description = "Real-space integration block size; zero selects it automatically.";
        item.default_value = "0";
        read_sync_int(input.bz);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const int values[3] = {para.input.bx, para.input.by, para.input.bz};
            for (int i = 0; i < 3; ++i)
            {
                if (values[i] < 0 || values[i] > 10)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "bx, by, and bz must be between 0 and 10.");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nb2d");
        item.annotation = "MPI 2D block size";
        item.category = "Parallel";
        item.type = "Integer";
        item.description = "Block size used to distribute H(R) in MPI builds; zero selects it automatically.";
        item.default_value = "0";
        read_sync_int(input.nb2d);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.nb2d < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nb2d must be non-negative.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nspin");
        item.annotation = "physical spin mode: 1, 2, or 4";
        item.category = "Density and spin";
        item.type = "Integer";
        item.description = "1 is spin-degenerate, 2 is collinear spin, and 4 is a complex two-component spinor.";
        item.default_value = "1";
        read_sync_int(input.nspin);
        item.reset_value = [](const Input_Item&, Parameter& para) {
            if (para.input.noncolin || para.input.lspinorb)
            {
                para.input.nspin = 4;
            }
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.nspin != 1 && para.input.nspin != 2 && para.input.nspin != 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nspin must be 1, 2, or 4.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("noncolin");
        item.annotation = "enable non-collinear magnetism";
        item.category = "Density and spin";
        item.type = "Boolean";
        item.description = "Enable non-collinear magnetism and select nspin=4.";
        item.default_value = "false";
        read_sync_bool(input.noncolin);
        this->add_item(item);
    }
    {
        Input_Item item("lspinorb");
        item.annotation = "enable spin-orbit coupling";
        item.category = "Density and spin";
        item.type = "Boolean";
        item.description = "Enable spin-orbit coupling, require a compatible fully relativistic pseudopotential, and select nspin=4.";
        item.default_value = "false";
        read_sync_bool(input.lspinorb);
        this->add_item(item);
    }
    {
        Input_Item item("soc_lambda");
        item.annotation = "SOC interpolation fraction";
        item.category = "Density and spin";
        item.type = "Real";
        item.description = "Scale between scalar-relativistic (0) and full pseudopotential SOC (1).";
        item.default_value = "1";
        read_sync_double(input.soc_lambda);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.soc_lambda < 0.0 || para.input.soc_lambda > 1.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "soc_lambda must be between 0 and 1.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("init_chg");
        item.annotation = "rho0 source for full H0";
        item.category = "Density and spin";
        item.type = "String";
        item.description = "atomic builds rho0 from isolated atoms; file requires existing restart/cube density; auto tries files first and falls back to atomic. Ignored by core.";
        item.default_value = "auto";
        read_sync_string(input.init_chg);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const std::vector<std::string> allowed = {"atomic", "file", "auto"};
            if (std::find(allowed.begin(), allowed.end(), para.input.init_chg) == allowed.end())
            {
                ModuleBase::WARNING_QUIT("ReadInput", "init_chg must be atomic, file, or auto.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dft_functional");
        item.annotation = "built-in LDA or GGA functional";
        item.category = "Density and spin";
        item.type = "String";
        item.description = "Use the pseudopotential functional by default. Explicit built-in choices are LDA/PZ, PWLDA, PBE, PBESOL, REVPBE, WC, BLYP, BP, PW91, HCTH, and OLYP. Libxc, meta-GGA, hybrid, HF, and EXX forms are unavailable.";
        item.default_value = "default";
        read_sync_string(input.dft_functional);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const std::vector<std::string> allowed = {"default", "lda", "pz", "slapznogxnogc", "pwlda",
                                                       "pbe", "slapwpbxpbc", "pbesol", "revpbe", "wc",
                                                       "blyp", "bp", "pw91", "hcth", "olyp"};
            if (std::find(allowed.begin(), allowed.end(), lowercase(para.input.dft_functional)) == allowed.end())
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dft_functional must select a built-in LDA or GGA functional.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nelec");
        item.annotation = "total electron count override";
        item.category = "Density and spin";
        item.type = "Real";
        item.description = "Override the electron count inferred from pseudopotentials; zero keeps the inferred value.";
        item.default_value = "0";
        read_sync_double(input.nelec);
        this->add_item(item);
    }
    {
        Input_Item item("nupdown");
        item.annotation = "spin-up minus spin-down electrons";
        item.category = "Density and spin";
        item.type = "Real";
        item.description = "Constrain the collinear initial spin imbalance used to construct rho0.";
        item.default_value = "0";
        read_sync_double(input.nupdown);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.nspin == 1 && para.input.nupdown != 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nupdown must be zero when nspin=1.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_rcut");
        item.annotation = "pseudopotential radial cutoff";
        item.category = "Density and spin";
        item.type = "Real";
        item.unit = "Bohr";
        item.description = "Maximum pseudopotential radial integration range.";
        item.default_value = "15";
        read_sync_double(input.pseudo_rcut);
        item.check_value = [](const Input_Item&, const Parameter& para) { require_positive("pseudo_rcut", para.input.pseudo_rcut); };
        this->add_item(item);
    }
    {
        Input_Item item("min_dist_coef");
        item.annotation = "minimum interatomic-distance factor";
        item.category = "Structure";
        item.type = "Real";
        item.description = "Reject structures with distances below this fraction of the reference covalent distance.";
        item.default_value = "0.2";
        read_sync_double(input.min_dist_coef);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.min_dist_coef < 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "min_dist_coef must be non-negative.");
            }
        };
        this->add_item(item);
    }
}

} // namespace ModuleIO
