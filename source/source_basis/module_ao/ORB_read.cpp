#include "ORB_read.h"

LCAO_Orbitals::LCAO_Orbitals() : Phi(new Numerical_Orbital[1]) {}

LCAO_Orbitals::~LCAO_Orbitals()
{
    delete[] Phi;
}

std::vector<double> LCAO_Orbitals::cutoffs() const
{
    std::vector<double> result(ntype);
    for (int it = 0; it < ntype; ++it)
    {
        result[it] = Phi[it].getRcut();
    }
    return result;
}
