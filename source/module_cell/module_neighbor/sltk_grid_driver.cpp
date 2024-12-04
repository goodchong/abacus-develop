#include "sltk_grid_driver.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"

#ifdef _OPENMP
#include <omp.h>
#endif

Grid_Driver::Grid_Driver(
	const int &test_d_in, 
	const int &test_grid_in)
:test_deconstructor(test_d_in),
Grid(test_grid_in)
{
	test_deconstructor	= test_d_in;
}

Grid_Driver::~Grid_Driver()
{
}

void Grid_Driver::Find_atom(const int& ntype,
                            const int& nnumber,
                            AdjacentAtomInfo* adjs) const
{
    ModuleBase::timer::tick("Grid_Driver", "Find_atom");
    //	std::cout << "lenght in Find atom = " << atomlink[offset].fatom.getAdjacentSet()->getLength() << std::endl;

    // store result in member adj_info when parameter adjs is NULL
    AdjacentAtomInfo* local_adjs = adjs == nullptr ? &this->adj_info : adjs;
    local_adjs->clear();
    const std::vector<FAtom*>& all_atom = all_adj_info[ntype][nnumber];

    for (const FAtom* atom: all_atom)
    {
        local_adjs->ntype.push_back(atom->type);
        local_adjs->natom.push_back(atom->natom);
        local_adjs->box.push_back(ModuleBase::Vector3<int>(atom->cell_x, atom->cell_y, atom->cell_z));
        local_adjs->adjacent_tau.push_back(ModuleBase::Vector3<double>(atom->x, atom->y, atom->z));
        local_adjs->adj_num++;
    }
    local_adjs->adj_num--;

    ModuleBase::timer::tick("Grid_Driver", "Find_atom");
    return;
}

// filter_adjs delete not adjacent atoms in adjs
void filter_adjs(const std::vector<bool>& is_adj, AdjacentAtomInfo& adjs)
{
    const int size = adjs.adj_num + 1;
    for (int i = size - 1; i >= 0; --i)
    {
        if (!is_adj[i])
        {
            adjs.adj_num--;
            adjs.ntype.erase(adjs.ntype.begin() + i);
            adjs.natom.erase(adjs.natom.begin() + i);
            adjs.adjacent_tau.erase(adjs.adjacent_tau.begin() + i); // info of adjacent_tau is not used in future
            adjs.box.erase(adjs.box.begin() + i);
        }
    }
}
