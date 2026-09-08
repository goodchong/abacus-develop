#ifndef UNITCELL_DATA_H
#define UNITCELL_DATA_H

#include <vector>

#include "source_base/intarray.h"
#include "source_base/matrix3.h"
/// @brief info of lattice
struct Lattice
{
    std::string Coordinate = "Direct"; ///< "Direct" or "Cartesian" or "Cartesian_angstrom"
    std::string latName = "user_defined_lattice";  ///< Lattice name
    double lat0 = 0.0;                 ///< Lattice constant(bohr)(a.u.)
    double lat0_angstrom = 0.0;        ///< Lattice constant(angstrom)
    double tpiba = 0.0;                ///< 2*pi / lat0;
    double tpiba2 = 0.0;               ///< tpiba ^ 2
    double omega = 0.0;                ///< the volume of the unit cell
    ModuleBase::Matrix3 latvec = ModuleBase::Matrix3();           ///< Unitcell lattice vectors
    ModuleBase::Vector3<double> a1, a2, a3;                       ///< Same as latvec, just at another form.
    ModuleBase::Vector3<double> latcenter;                        ///< (a1+a2+a3)/2 the center of vector
    ModuleBase::Matrix3 G = ModuleBase::Matrix3();                ///< reciprocal lattice vector (2pi*inv(R) )
    ModuleBase::Matrix3 GT = ModuleBase::Matrix3();               ///< traspose of G
    ModuleBase::Matrix3 GGT = ModuleBase::Matrix3();              ///< GGT = G*GT
    ModuleBase::Matrix3 invGGT = ModuleBase::Matrix3();           ///< inverse G
};

//========================================================
// relationship between:
// ntype, it
// nat, iat
// atoms[it].na, ia,
// atoms[it].nw, iw
//
// if know it ==> atoms[it].na; atoms[it].nw
// if know iat ==> it; ia;
// if know ia, mush have known it ==> iat
// if know iwt, must have known it, ia ==> iwt
//========================================================
/// @brief usefull data and index maps
struct Statistics
{
    int ntype = 0;                 // number of atom species in UnitCell
    int nat = 0;                   // total number of atoms of all species in unitcell
    int* iat2it = nullptr;         // iat==>it, distinguish a atom belong to which type
    int* iat2ia = nullptr;         // iat==>ia
    ModuleBase::IntArray itia2iat; //(it, ia)==>iat, the index in nat, add 2009-3-2 by mohan
    int namax = 0;                 // the max na among all atom species

    ~Statistics()
    {
        delete[] iat2it;
        delete[] iat2ia;
    }
};

#endif
