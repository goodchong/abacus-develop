#ifndef GRID_H
#define GRID_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"
#include "sltk_util.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <list>

class Grid
{
  public:
    // Constructors and destructor
    // Grid is Global class,so init it with constant number
    Grid() : test_grid(0){};
    Grid(const int& test_grid_in);
    virtual ~Grid();

    void init(std::ofstream& ofs, const UnitCell& ucell, const double radius_in, const bool boundary = true);

    // Data
    bool pbc; // When pbc is set to false, periodic boundary conditions are explicitly ignored.
    double sradius2; // searching radius squared (unit:lat0)
    double sradius;  // searching radius (unit:lat0)
    
    // coordinate range of the input atom (unit:lat0)
    double x_min;
    double y_min;
    double z_min;
    double x_max;
    double y_max;
    double z_max;

    // The algorithm for searching neighboring atoms uses a "box" partitioning method. 
    // Each box has an edge length of sradius, and the number of boxes in each direction is recorded here.
    double box_edge_length;
    int box_nx;
    int box_ny;
    int box_nz;

    void getBox(int& bx, int& by, int& bz, const double& x, const double& y, const double& z)
    {
        bx = std::floor((x - x_min) / box_edge_length);
        by = std::floor((y - y_min) / box_edge_length);
        bz = std::floor((z - z_min) / box_edge_length);
    }

    // 20241210 zhanghaochong 
    // Here, instances of all atoms are stored as a single copy in the atoms_in_box structure, 
    // which is a quadruple nested combination of vector and list. The all_adj_info structure only 
    // stores pointers that reference the memory of atoms_in_box. Please note that vector is !NOT! a 
    // memory-stable data structure. During operations such as resize or push_back, the addresses of 
    // existing data within the vector may change. Therefore, if atoms_in_box undergoes resize, push_back, 
    // or erase after all_adj_info has been set to point to it, the pointers in all_adj_info will 
    // result in memory leaks.
    // Therefore, atoms_in_box and all_adj_info vector part must not undergo any push_back or erase operations during 
    // use. After resizing during initialization, they should not be modified.

    // Stores the atoms after box partitioning.
    std::vector<std::vector<std::vector< std::list<FAtom> >>> atoms_in_box;

    // Stores the adjacent information of atoms. [ntype][natom][adj list]
    std::vector<std::vector< std::list<FAtom *> >> all_adj_info;
    // There are two main benefits to storing atomic information in a single copy and using pointers for 
    // the other. On one hand, the construction and copying cost of the fatom class during computation 
    // should not be underestimated. On the other hand, it saves memory. If we need to deal with millions 
    // of atoms, just storing the atoms themselves would consume a significant amount of memory.

    void clear_atoms()
    {
        // we have to clear the all_adj_info
        // because the pointers point to the memory in vector atoms_in_box
        all_adj_info.clear();
        atoms_in_box.clear();
    }
    void clear_adj_info()
    {
        // here dont need to free the memory, 
        // because the pointers point to the memory in vector atoms_in_box
        all_adj_info.clear();
    }
    int getGlayerX() const
    {
        return glayerX;
    }
    int getGlayerY() const
    {
        return glayerY;
    }
    int getGlayerZ() const
    {
        return glayerZ;
    }
    int getGlayerX_minus() const
    {
        return glayerX_minus;
    }
    int getGlayerY_minus() const
    {
        return glayerY_minus;
    }
    int getGlayerZ_minus() const
    {
        return glayerZ_minus;
    }
  private:
    const int test_grid;

    void setMemberVariables(std::ofstream& ofs_in, const UnitCell& ucell);

    void Construct_Adjacent(const UnitCell& ucell);
    void Construct_Adjacent_near_box(const FAtom& fatom);
    void Construct_Adjacent_final(const FAtom& fatom1, FAtom* fatom2);

    void Check_Expand_Condition(const UnitCell& ucell);
    int glayerX;
    int glayerX_minus;
    int glayerY;
    int glayerY_minus;
    int glayerZ;
    int glayerZ_minus;
};

#endif
