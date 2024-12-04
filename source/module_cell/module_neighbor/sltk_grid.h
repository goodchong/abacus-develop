#ifndef GRID_H
#define GRID_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"
#include "sltk_atom_input.h"
#include "sltk_util.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

typedef std::vector<FAtom> AtomMap;

class Grid
{
  public:
    // Constructors and destructor
    // Grid is Global class,so init it with constant number
    Grid() : test_grid(0){};
    Grid(const int& test_grid_in);
    virtual ~Grid();

    void init(std::ofstream& ofs, const UnitCell& ucell, const Atom_input& input);

    // 2015-05-07
    void delete_vector(int i, int j, int k);

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

    // If there is no cells expansion, this would be 1, 1, 1.
    // If a cell is expanded, it indicates the number of unit cells in each direction, 
    // including the original unit cell. For example, 3, 3, 3 would mean 27 unit cells.
    int cell_nx;
    int cell_ny;
    int cell_nz;

    // true cell means the index of original cell. cell index start from 0 to (nx-1)
    int true_cell_x;
    int true_cell_y;
    int true_cell_z;

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
    // Stores the atoms after box partitioning.
    std::vector<std::vector<std::vector<AtomMap>>> atoms_in_box;

    // Stores the adjacent information of atoms. [ntype][natom][adj list]
    std::vector<std::vector< std::vector<FAtom *> >> all_adj_info;
    void clear_atoms()
    {
        atoms_in_box.clear();
    }

    // LiuXh add 2019-07-15
    int getCellX() const
    {
        return cell_nx;
    }
    int getCellY() const
    {
        return cell_ny;
    }
    int getCellZ() const
    {
        return cell_nz;
    }
    int getTrueCellX() const
    {
        return true_cell_x;
    }
    int getTrueCellY() const
    {
        return true_cell_y;
    }
    int getTrueCellZ() const
    {
        return true_cell_z;
    }

  private:
    const int test_grid;

    void setMemberVariables(std::ofstream& ofs_in, const UnitCell& ucell, const Atom_input& input);

    void setBoundaryAdjacent(std::ofstream& ofs_in, const Atom_input& input);

    void Construct_Adjacent(const UnitCell& ucell);

    void Construct_Adjacent_expand_periodic(FAtom& fatom);

    void Construct_Adjacent_final(FAtom& fatom1,
                                  FAtom& fatom2);
};

#endif
