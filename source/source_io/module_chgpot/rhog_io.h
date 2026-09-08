#ifndef RHOG_IO_H
#define RHOG_IO_H

#include <complex>
#include <string>

#include "source_basis/module_pw/pw_basis.h"

/**
 * Read the legacy ABACUS binary rho(G) restart format. The stored header
 * includes the reciprocal-grid representation flag, G-vector count, spin
 * count, lattice vectors, and Miller indices before the complex rho(G) data.
 */

namespace ModuleIO
{

bool read_rhog(const std::string& filename, const ModulePW::PW_Basis* pw_rhod, std::complex<double>** rhog);

} // namespace ModuleIO

#endif
