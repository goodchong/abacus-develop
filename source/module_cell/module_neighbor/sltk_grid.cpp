#include "sltk_grid.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "sltk_atom_input.h"

Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
{
}

Grid::~Grid()
{
    this->clear_atoms();
}

void Grid::init(std::ofstream& ofs_in, const UnitCell& ucell, const Atom_input& input)
{
    ModuleBase::TITLE("SLTK_Grid", "init");
    ModuleBase::timer::tick("atom_arrange", "grid_d.init");

    this->setMemberVariables(ofs_in, ucell, input);
    this->Construct_Adjacent(ucell);
    ModuleBase::timer::tick("atom_arrange", "grid_d.init");
}
void Grid::setMemberVariables(std::ofstream& ofs_in, //  output data to ofs
                              const UnitCell& ucell,
                              const Atom_input& input)
{
    ModuleBase::TITLE("SLTK_Grid", "setMemberVariables");

    this->clear_atoms();

    this->pbc = input.getBoundary();
    this->sradius2 = input.getRadius() * input.getRadius();
    this->sradius = input.getRadius();

    if (test_grid)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "PeriodicBoundary", this->pbc);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Radius(unit:lat0)", sradius);
    }

    // random selection, in order to estimate again.
    this->x_min = ucell.atoms[0].tau[0].x;
    this->y_min = ucell.atoms[0].tau[0].y;
    this->z_min = ucell.atoms[0].tau[0].z;
    this->x_max = ucell.atoms[0].tau[0].x;
    this->y_max = ucell.atoms[0].tau[0].y;
    this->z_max = ucell.atoms[0].tau[0].z;

    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    // calculate min & max value
    for (int ix = -input.getGrid_layerX_minus(); ix < input.getGrid_layerX(); ix++)
    {
        for (int iy = -input.getGrid_layerY_minus(); iy < input.getGrid_layerY(); iy++)
        {
            for (int iz = -input.getGrid_layerZ_minus(); iz < input.getGrid_layerZ(); iz++)
            {
                for (int i = 0; i < ucell.ntype; i++)
                {
                    for (int j = 0; j < ucell.atoms[i].na; j++)
                    {
                        double x = ucell.atoms[i].tau[j].x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double y = ucell.atoms[i].tau[j].y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double z = ucell.atoms[i].tau[j].z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
                        x_min = std::min(x_min, x);
                        x_max = std::max(x_max, x);
                        y_min = std::min(y_min, y);
                        y_max = std::max(y_max, y);
                        z_min = std::min(z_min, z);
                        z_max = std::max(z_max, z);
                    }
                }
            }
        }
    }
    ModuleBase::GlobalFunc::OUT(ofs_in, "Find the coordinate range of the input atom(unit:lat0).");
    ModuleBase::GlobalFunc::OUT(ofs_in, "min_tau", x_min, y_min, z_min);
    ModuleBase::GlobalFunc::OUT(ofs_in, "max_tau", x_max, y_max, z_max);

    this->cell_nx = input.getGrid_layerX() + input.getGrid_layerX_minus();
    this->cell_ny = input.getGrid_layerY() + input.getGrid_layerY_minus();
    this->cell_nz = input.getGrid_layerZ() + input.getGrid_layerZ_minus();
    this->true_cell_x = input.getGrid_layerX_minus();
    this->true_cell_y = input.getGrid_layerY_minus();
    this->true_cell_z = input.getGrid_layerZ_minus();

    this->box_edge_length = sradius + 0.1; // To avoid edge cases, the size of the box is slightly increased.

    this->box_nx = std::ceil((this->x_max - this->x_min) / box_edge_length) + 1;
    this->box_ny = std::ceil((this->y_max - this->y_min) / box_edge_length) + 1;
    this->box_nz = std::ceil((this->z_max - this->z_min) / box_edge_length) + 1;

    ModuleBase::GlobalFunc::OUT(ofs_in, "CellNumber", cell_nx, cell_ny, cell_nz);
    ModuleBase::GlobalFunc::OUT(ofs_in, "BoxNumber", box_nx, box_ny, box_nz);
    ModuleBase::GlobalFunc::OUT(ofs_in, "TrueCellNumber", true_cell_x, true_cell_y, true_cell_z);

    atoms_in_box.resize(this->box_nx);
    for (int i = 0; i < this->box_nx; i++)
    {
        atoms_in_box[i].resize(this->box_ny);
        for (int j = 0; j < this->box_ny; j++)
        {
            atoms_in_box[i][j].resize(this->box_nz);
        }
    }


    for (int ix = -input.getGrid_layerX_minus(); ix < input.getGrid_layerX(); ix++)
    {
        for (int iy = -input.getGrid_layerY_minus(); iy < input.getGrid_layerY(); iy++)
        {
            for (int iz = -input.getGrid_layerZ_minus(); iz < input.getGrid_layerZ(); iz++)
            {
                for (int i = 0; i < ucell.ntype; i++)
                {
                    for (int j = 0; j < ucell.atoms[i].na; j++)
                    {
                        double x = ucell.atoms[i].tau[j].x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double y = ucell.atoms[i].tau[j].y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double z = ucell.atoms[i].tau[j].z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
                        FAtom atom(x, y, z, i, j, ix, iy, iz);
                        int box_i_x, box_i_y, box_i_z;
                        this->getBox(box_i_x, box_i_y, box_i_z, x, y, z);
                        this->atoms_in_box[box_i_x][box_i_y][box_i_z].push_back(atom);
                    }
                }
            }
        }
    }
    
    this->all_adj_info.resize(ucell.ntype);
    for (int i = 0; i < ucell.ntype; i++)
    {
        this->all_adj_info[i].resize(ucell.atoms[i].na);
    }
}

void Grid::Construct_Adjacent(const UnitCell& ucell)
{
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");

    for  (int i_type = 0; i_type < ucell.ntype; i_type++)
    {
        for (int j_atom = 0; j_atom < ucell.atoms[i_type].na; j_atom++)
        {

            FAtom atom(ucell.atoms[i_type].tau[j_atom].x,
                     ucell.atoms[i_type].tau[j_atom].y,
                     ucell.atoms[i_type].tau[j_atom].z,
                     i_type,
                     j_atom,
                     0, 0 ,0);

            this->Construct_Adjacent_expand_periodic(atom);
        }
    }
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");
}

void Grid::Construct_Adjacent_expand_periodic(FAtom& fatom)
{
    //	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_expand_periodic");
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");
    int box_i_x, box_i_y, box_i_z;
    this->getBox(box_i_x, box_i_y, box_i_z, fatom.x(), fatom.y(), fatom.z());

    for (int box_i_x_adj = std::max(box_i_x - 1, 0); box_i_x_adj <= std::min(box_i_x + 1, box_nx - 1); box_i_x_adj++)
    {
        for (int box_i_y_adj = std::max(box_i_y - 1, 0); box_i_y_adj <= std::min(box_i_y + 1, box_ny - 1); box_i_y_adj++)
        {
            for (int box_i_z_adj = std::max(box_i_z - 1, 0); box_i_z_adj <= std::min(box_i_z + 1, box_nz - 1); box_i_z_adj++)
            {
                for (auto &fatom2 : this->atoms_in_box[box_i_x_adj][box_i_y_adj][box_i_z_adj])
                {
                    this->Construct_Adjacent_final(fatom, fatom2);
                }
            }
        }
    }
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");
}

void Grid::Construct_Adjacent_final(FAtom& fatom1,
                                    FAtom& fatom2)
{
    double delta_x = fatom1.x() - fatom2.x();
    double delta_y = fatom1.y() - fatom2.y();
    double delta_z = fatom1.z() - fatom2.z();

    double dr = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

    if (dr <= this->sradius2)
    {
        all_adj_info[fatom1.getType()][fatom1.getNatom()].push_back(&fatom2);
    }
}
