#ifndef WRITE_H0_H
#define WRITE_H0_H

#include <string>

class Parallel_Orbitals;
class UnitCell;

namespace hamilt
{
template <typename T>
class HContainer;
}

namespace ModuleIO
{

/**
 * @brief Gather and write one spin channel of the H0 real-space matrix.
 *
 * The numeric header layout and CSR body match the established ABACUS text
 * H(R) format.  H0-specific metadata is confined to the existing comment
 * line so readers which consume the numeric fields are unaffected.
 */
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
              const int rank);

} // namespace ModuleIO

#endif // WRITE_H0_H
