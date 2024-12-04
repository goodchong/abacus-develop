#include "sltk_atom_input.h"

#include "module_base/memory.h"
#include "module_parameter/parameter.h"
#include "sltk_grid.h"

//==========================================================
// define constructor and deconstructor
//==========================================================
Atom_input::Atom_input(std::ofstream& ofs_in,
                       const UnitCell& ucell,
                       const bool boundary_in,
                       const double radius_in,
                       const int& test_atom_in)
    : periodic_boundary(boundary_in), radius(radius_in),
      glayerX(1), glayerX_minus(0), glayerY(1),
      glayerY_minus(0), glayerZ(1), glayerZ_minus(0),
      test_atom_input(test_atom_in)
{
    ModuleBase::TITLE("Atom_input", "Atom_input");

    if (test_atom_input)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "Periodic_boundary", periodic_boundary);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Searching radius(lat0)", radius);
    }

    this->Check_Expand_Condition(ucell);

    if (test_atom_input)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "glayer+", glayerX, glayerY, glayerZ);
        ModuleBase::GlobalFunc::OUT(ofs_in, "glayer-", glayerX_minus, glayerY_minus, glayerZ_minus);
    }
}

Atom_input::~Atom_input()
{
}

//============================================
// !!!! May still have bug, be very careful!!
// should use the same algorithm to generate
// dxe, dye, dze in grid_meshcell.cpp.
//============================================
void Atom_input::Check_Expand_Condition(const UnitCell& ucell)
{
    //	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

    if (!periodic_boundary)
    {
        return;
    }

    /*2016-07-19, LiuXh
        // the unit of extent_1DX,Y,Z is lat0.
        // means still how far can be included now.
        double extent_1DX = glayerX * clength0 - dmaxX;
        while (radius > extent_1DX)
        {
            glayerX++;
            extent_1DX = glayerX * clength0 - dmaxX;
        }
        double extent_1DY = glayerY * clength1 - dmaxY;
        while (radius > extent_1DY)
        {
            glayerY++;
            extent_1DY = glayerY * clength1 - dmaxY;
        }
        double extent_1DZ = glayerZ * clength2 - dmaxZ;
        while (radius > extent_1DZ)
        {
            glayerZ++;
            extent_1DZ = glayerZ * clength2 - dmaxZ;
        }

        // in case the cell is not retangle.
        // mohan added 2009-10-23
        // if this is not added, it's a serious bug.
        glayerX++;
        glayerY++;
        glayerZ++;
        if(test_atom_input)
        {
            GlobalV::ofs_running << " Extend distance from the (maxX,maxY,maxZ) direct position in this unitcell: " <<
    std::endl;
        }

        if(test_atom_input)OUT(GlobalV::ofs_running,"ExtentDim+",extent_1DX,extent_1DY,extent_1DZ);

        double extent_1DX_minus = glayerX_minus * clength0 + dminX;
        while (radius > extent_1DX_minus)
        {
            glayerX_minus++;
            extent_1DX_minus = glayerX_minus * clength0 + dminX;
        }
        double extent_1DY_minus = glayerY_minus * clength1 + dminY;
        while (radius > extent_1DY_minus)
        {
            glayerY_minus++;
            extent_1DY_minus = glayerY_minus * clength1 + dminY;
        }
        double extent_1DZ_minus = glayerZ_minus * clength2 + dminZ;
        while (radius > extent_1DZ_minus)
        {
            glayerZ_minus++;
            extent_1DZ_minus = glayerZ_minus * clength2 + dminZ;
        }

        // in case the cell is not retangle.
        // mohan added 2009-10-23
        // if this is not added, it's a serious bug.
        glayerX_minus++;
        glayerY_minus++;
        glayerZ_minus++;

        //glayerX_minus++;
        //glayerY_minus++;
        //glayerZ_minus++;
    2016-07-19, LiuXh*/
    // Begin, 2016-07-19, LiuXh
    double a23_1 = ucell.latvec.e22 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e32;
    double a23_2 = ucell.latvec.e21 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e31;
    double a23_3 = ucell.latvec.e21 * ucell.latvec.e32 - ucell.latvec.e22 * ucell.latvec.e31;
    double a23_norm = sqrt(a23_1 * a23_1 + a23_2 * a23_2 + a23_3 * a23_3);
    double extend_v = a23_norm * radius;
    double extend_d1 = extend_v / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d11 = std::ceil(extend_d1);

    double a31_1 = ucell.latvec.e32 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e12;
    double a31_2 = ucell.latvec.e31 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e11;
    double a31_3 = ucell.latvec.e31 * ucell.latvec.e12 - ucell.latvec.e32 * ucell.latvec.e11;
    double a31_norm = sqrt(a31_1 * a31_1 + a31_2 * a31_2 + a31_3 * a31_3);
    double extend_d2 = a31_norm * radius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d22 = std::ceil(extend_d2);

    double a12_1 = ucell.latvec.e12 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e22;
    double a12_2 = ucell.latvec.e11 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e21;
    double a12_3 = ucell.latvec.e11 * ucell.latvec.e22 - ucell.latvec.e12 * ucell.latvec.e21;
    double a12_norm = sqrt(a12_1 * a12_1 + a12_2 * a12_2 + a12_3 * a12_3);
    double extend_d3 = a12_norm * radius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d33 = std::ceil(extend_d3);
    // 2016-09-05, LiuXh

    glayerX = extend_d11 + 1;
    glayerY = extend_d22 + 1;
    glayerZ = extend_d33 + 1;
    glayerX_minus = extend_d11;
    glayerY_minus = extend_d22;
    glayerZ_minus = extend_d33;
    // End, 2016-09-05, LiuXh
}
