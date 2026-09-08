#include "write_h0.h"

#include "source_base/tool_quit.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_io/module_output/ucell_io.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"

#include <complex>
#include <fstream>

namespace ModuleIO
{

namespace
{

template <typename T>
void write_h0_serial(hamilt::HContainer<T>* h0,
                     const UnitCell& ucell,
                     const std::string& filename,
                     const std::string& h0_type,
                     const int physical_nspin,
                     const int output_nspin,
                     const int ispin,
                     const double sparse_threshold,
                     const int precision)
{
    std::ofstream ofs(filename.c_str());
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_h0", "Cannot open H0 CSR file: " + filename);
    }

    const char* formula = h0_type == "core"
                              ? "T + V_nl + V_loc"
                              : "T + V_nl + V_loc + V_H[rho0] + V_xc[rho0]";
    ofs << " --- Ionic Step 1 ---\n";
    ofs << " # print H0 matrix in real space H0(R); h0_type=" << h0_type
        << "; formula=" << formula << "; physical_nspin=" << physical_nspin
        << "; unit=Ry\n";
    ofs << " " << output_nspin << " # number of spin directions\n";
    ofs << " " << ispin + 1 << " # spin index\n";
    ofs << " " << h0->get_nbasis() << " # number of localized basis\n";
    ofs << " " << h0->size_R_loop() << " # number of Bravais lattice vector R\n\n";

    UcellIO::write_ucell(ofs, &ucell);
    ofs << '\n';

    hamilt::Output_HContainer<T> output(h0, ofs, sparse_threshold, precision);
    output.write();
}

} // namespace

template <typename T>
void write_h0(const hamilt::HContainer<T>& h0,
              const UnitCell& ucell,
              const Parallel_Orbitals& para_v,
              const std::string& output_dir,
              const std::string& h0_type,
              const int physical_nspin,
              const int output_nspin,
              const int ispin,
              const double sparse_threshold,
              const int precision,
              const int rank)
{
    const int nbasis = h0.get_nbasis();
#ifdef __MPI
    Parallel_Orbitals serial_v;
    serial_v.init(nbasis, nbasis, nbasis, para_v.comm());
    serial_v.set_serial(nbasis, nbasis);
    serial_v.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, nbasis);
    hamilt::HContainer<T> h0_serial(&serial_v);
    hamilt::gatherParallels(h0, &h0_serial, 0);
#else
    hamilt::HContainer<T> h0_serial(h0);
    h0_serial.add_value_intersection(h0);
#endif

    if (rank == 0)
    {
        const std::string filename = output_dir + "hrs" + std::to_string(ispin + 1) + "_nao.csr";
        write_h0_serial(&h0_serial,
                        ucell,
                        filename,
                        h0_type,
                        physical_nspin,
                        output_nspin,
                        ispin,
                        sparse_threshold,
                        precision);
    }
}

template void write_h0<double>(const hamilt::HContainer<double>&,
                               const UnitCell&,
                               const Parallel_Orbitals&,
                               const std::string&,
                               const std::string&,
                               const int,
                               const int,
                               const int,
                               const double,
                               const int,
                               const int);
template void write_h0<std::complex<double>>(const hamilt::HContainer<std::complex<double>>&,
                                             const UnitCell&,
                                             const Parallel_Orbitals&,
                                             const std::string&,
                                             const std::string&,
                                             const int,
                                             const int,
                                             const int,
                                             const double,
                                             const int,
                                             const int);

} // namespace ModuleIO
