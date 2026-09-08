#include <cstring>        // Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include <regex>
#include <fstream>

#include "unitcell.h"
#include "read_atoms_helper.h"

#include "print_cell.h"
#include "read_stru.h"
#include "source_base/timer.h"
#include "source_base/constants.h"
#include "source_base/formatter.h"
#include "source_base/mathzone.h"

bool unitcell::read_atom_positions(UnitCell& ucell,
                         std::ifstream &ifpos,
                         std::ofstream &ofs_running,
                         std::ofstream &ofs_warning,
                         const int nspin,
                         const std::string& orbital_dir,
                         const bool noncolin)
{
    ModuleBase::TITLE("UnitCell","read_atom_positions");

    std::string& Coordinate  = ucell.Coordinate;
    const int ntype = ucell.ntype;
    assert (nspin==1 || nspin==2 || nspin==4);

    if (ucell.magnet.start_mag.size() != static_cast<size_t>(ntype))
    {
        ucell.magnet.start_mag.resize(ntype, 0.0);
    }

    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifpos, "ATOMIC_POSITIONS"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifpos, Coordinate);

        if (!unitcell::validate_coordinate_system(Coordinate, ofs_warning))
        {
            return false;
        }

        ucell.nat = 0;

        //======================================
        // calculate total number of ucell.atoms
        // and adjust the order of atom species
        //======================================
        for (int it = 0;it < ntype; it++)
        {
            ofs_running << "\n READING ATOM TYPE " << it+1 << std::endl;

            bool set_element_mag_zero = false;
            if (!unitcell::read_atom_type_header(it, ucell, ifpos, ofs_running,
                                       ofs_warning, set_element_mag_zero,
                                       orbital_dir))
            {
                return false;
            }

            int na = ucell.atoms[it].na;
            ucell.nat += na;

            if (na > 0)
            {
                unitcell::allocate_atom_properties(ucell.atoms[it], na, ucell.atom_mass[it]);
                for (int ia = 0;ia < na; ia++)
                {
                    ModuleBase::Vector3<double> v;
                    ifpos >> v.x >> v.y >> v.z;
                    ucell.atoms[it].mag[ia]=ucell.magnet.start_mag[it];
                    //if this line is used, default startmag_type would be 2
                    ucell.atoms[it].angle1[ia]=0;
                    ucell.atoms[it].angle2[ia]=0;
                    ucell.atoms[it].m_loc_[ia].set(0,0,0);

                    bool input_vec_mag=false;
                    bool input_angle_mag=false;

                    // Parse optional properties
                    if (!unitcell::parse_atom_properties(ifpos, ucell.atoms[it], ia,
                                              input_vec_mag, input_angle_mag,
                                              set_element_mag_zero))
                    {
                        return false;
                    }

                    // Process magnetization
                    unitcell::process_magnetization(ucell.atoms[it], it, ia, nspin,
                                        input_vec_mag, input_angle_mag, ofs_running,
                                        noncolin);

                    // Transform coordinates
                    unitcell::transform_atom_coordinates(ucell.atoms[it], ia, Coordinate,
                                             v, ucell.latvec, ucell.lat0, ucell.latcenter);

                }//endj
            }    // end na
            // reset some useless parameters
            if (set_element_mag_zero)
            {
                ucell.magnet.start_mag[it] = 0.0;
            }
        } // end for ntype

        // Auto-set magnetization if needed
        unitcell::autoset_magnetization(ucell, nspin, ofs_running);
    }   // end scan_begin

    // Final validation and output
    return unitcell::finalize_atom_positions(ucell, ofs_running, ofs_warning);

}//end read_atom_positions
