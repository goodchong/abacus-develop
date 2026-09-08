#ifndef UCELL_IO_H
#define UCELL_IO_H

#include "source_cell/unitcell.h"

#include <fstream>

namespace ModuleIO {

/**
 * @brief A class for unit cell I/O operations
 * 
 * This class writes the legacy cell preamble retained by the H(R) CSR format.
 */
class UcellIO {
public:
    /**
     * @brief Writes the unit cell information to a file.
     *
     * @param ofs The output file stream.
     * @param ucell A pointer to the UnitCell object.
     */
    static void write_ucell(std::ofstream& ofs, const UnitCell* ucell);

};

} // namespace ModuleIO

#endif // UCELL_IO_H
