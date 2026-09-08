#include <regex>
#include <cassert>

#include "print_cell.h"
#include "source_base/formatter.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"

namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0,
                   std::ofstream &ofs)
    {
        ModuleBase::TITLE("UnitCell", "print_tau");
        // assert (direct || Coordinate == "Cartesian" || Coordinate == "Cartesian_angstrom"); // this line causes abort in unittest ReadAtomPositionsCACXY.
        // previously there are two if-statements, the first is `if(Coordinate == "Direct")` and the second is `if(Coordinate == "Cartesian" || Coordiante == "Cartesian_angstrom")`
        // however the Coordinate can also be value among Cartesian_angstrom_center_xy, Cartesian_angstrom_center_xz, Cartesian_angstrom_center_yz and Cartesian_angstrom_center_xyz

        // if Coordinate has value one of them, this print_tau will not print anything.
        std::regex pattern("Direct|Cartesian(_angstrom)?(_center_(xy|xz|yz|xyz))?");
        assert(std::regex_search(Coordinate, pattern));
        bool direct = (Coordinate == "Direct");

        //----------------------
        // print atom positions
        //----------------------
        std::string table;
        table += direct? " DIRECT COORDINATES\n": FmtCore::format(" CARTESIAN COORDINATES ( UNIT = %15.8f Bohr )\n", lat0);
        table += FmtCore::format("%5s%19s%19s%19s%8s\n", "atom", "x", "y", "z", "mag");
        for(int it = 0; it < ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
                const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
                const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
                table += FmtCore::format("%5s%19.12f%19.12f%19.12f%8.4f\n", 
                                        atoms[it].label, 
                                        x, 
                                        y, 
                                        z, 
                                        atoms[it].mag[ia]); 
            }
        }
        table += "\n";
        ofs << table; 

        return;
    }

}
