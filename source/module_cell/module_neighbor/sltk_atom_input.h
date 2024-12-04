#ifndef ATOM_INPUT_H
#define ATOM_INPUT_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"

class Atom_input
{
  public:
    //==========================================================
    // Constructors and destructor
    //==========================================================
    Atom_input(std::ofstream& ofs_in,
               const UnitCell& ucell,
               const bool boundary = true, // 1 : periodic ocndition
               const double radius_in = 0, // searching radius
               const int& test_atom_in = 0 // caoyu reconst 2021-05-24
    );
    ~Atom_input();

  public:
    int getBoundary() const
    {
        return periodic_boundary;
    }

    double getRadius() const
    {
        return radius;
    }

    int getGrid_layerX() const
    {
        return glayerX;
    }

    int getGrid_layerX_minus() const
    {
        return glayerX_minus;
    }

    int getGrid_layerY() const
    {
        return glayerY;
    }

    int getGrid_layerY_minus() const
    {
        return glayerY_minus;
    }

    int getGrid_layerZ() const
    {
        return glayerZ;
    }

    int getGrid_layerZ_minus() const
    {
        return glayerZ_minus;
    }

  private:
    int test_atom_input; // caoyu reconst 2021-05-24
    bool periodic_boundary;
    double radius;

    void Check_Expand_Condition(const UnitCell& ucell);
    int glayerX;
    int glayerX_minus;
    int glayerY;
    int glayerY_minus;
    int glayerZ;
    int glayerZ_minus;
};

#endif
