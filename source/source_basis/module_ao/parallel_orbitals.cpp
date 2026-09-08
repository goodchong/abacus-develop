#include "parallel_orbitals.h"

#include "source_base/global_function.h"

Parallel_Orbitals::Parallel_Orbitals()
{
    this->nnr = 1;
}

Parallel_Orbitals::~Parallel_Orbitals()
{
}

void Parallel_Orbitals::set_atomic_trace(const int* iat2iwt, const int &nat, const int &nlocal)
{
    ModuleBase::TITLE("Parallel_Orbitals", "set_atomic_trace");
    this->iat2iwt_ = iat2iwt;
    int nat_plus_1 = nat + 1;
    this->atom_begin_col.resize(nat_plus_1);
    this->atom_begin_row.resize(nat_plus_1);
    for(int iat=0;iat<nat;iat++)
    {
        this->atom_begin_col[iat] = -1;
        this->atom_begin_row[iat] = -1;
        int irow = iat2iwt[iat];
        int icol = iat2iwt[iat];
        const int nw_global = (iat == nat-1) ? (nlocal - irow): (iat2iwt[iat+1] - irow);
        //find the first local row index of atom iat
        for(int i=0;i<nw_global;i++)
        {
            if (this->global2local_row_[irow] != -1)
            {
                this->atom_begin_row[iat] = this->global2local_row_[irow];
                break;
            }
            irow++;
        }
        //find the first local col index of atom iat
        for(int i=0;i<nw_global;i++)
        {
            if (this->global2local_col_[icol] != -1)
            {
                this->atom_begin_col[iat] = this->global2local_col_[icol];
                break;
            }
            icol++;
        }
    }
    this->atom_begin_row[nat] = this->nrow;
    this->atom_begin_col[nat] = this->ncol;
}

// Get the number of columns of the orbital matrix of the iat-th atom
int Parallel_Orbitals::get_ncol_atom(int iat) const
{
    int size = this->atom_begin_col[iat];
    // If the iat-th atom does not have an orbital matrix, return 0
    if(size == -1)
    {
        return 0;
    }
    iat += 1;
    // Traverse the orbital matrices of the atom and calculate the number of columns
    while(this->atom_begin_col[iat] <= this->ncol)
    {
        if(this->atom_begin_col[iat] != -1)
        {
            size = this->atom_begin_col[iat] - size;
            return size;
        }
        iat++;
    }
    // If the orbital matrix is not found after all atoms are traversed, throw an exception
    throw std::string("error in get_ncol_atom(iat)");
}
// Get the number of rows of the orbital matrix of the iat-th atom
int Parallel_Orbitals::get_nrow_atom(int iat) const
{
    int size = this->atom_begin_row[iat];
    if(size == -1)
    {
        return 0;
    }
    iat += 1;
    while(this->atom_begin_row[iat] <= this->nrow)
    {
        if(this->atom_begin_row[iat] != -1)
        {
            size = this->atom_begin_row[iat] - size;
            return size;
        }
        iat++;
    }
    // If the orbital matrix is not found after all atoms are traversed, throw an exception
    throw std::string("error in get_nrow_atom(iat)");
}

bool Parallel_Orbitals::is_invalid_atom_pair(int iat1, int iat2) const
{
    return get_nrow_atom(iat1) <= 0 || get_ncol_atom(iat2) <= 0;
}

// Get the global indexes of the rows of the parallel orbital matrix
std::vector<int> Parallel_Orbitals::get_indexes_row() const
{
    std::vector<int> indexes(this->nrow);
    for(int i = 0; i < this->nrow; i++)
    {
#ifdef __MPI
        indexes[i] = this->local2global_row(i);
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
// Get the global indexes of the columns of the parallel orbital matrix
std::vector<int> Parallel_Orbitals::get_indexes_col() const
{
    std::vector<int> indexes(this->ncol);
    for(int i = 0; i < this->ncol; i++)
    {
#ifdef __MPI
        indexes[i] = this->local2global_col(i);
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
// Get the global indexes of the rows of the orbital matrix of the iat-th atom
std::vector<int> Parallel_Orbitals::get_indexes_row(int iat) const
{
    int size = this->get_nrow_atom(iat);
    if(size == 0)
    {
        return std::vector<int>();
    }
    std::vector<int> indexes(size);
    int irow = this->atom_begin_row[iat];
    int begin = this->iat2iwt_[iat];
    for(int i = 0; i < size; ++i)
    {
#ifdef __MPI
        indexes[i] = this->local2global_row(irow + i) - begin;
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
// Get the global indexes of the columns of the orbital matrix of the iat-th atom
std::vector<int> Parallel_Orbitals::get_indexes_col(int iat) const
{
    int size = this->get_ncol_atom(iat);
    if(size == 0)
    {
        return std::vector<int>();
    }
    std::vector<int> indexes(size);
    int icol = this->atom_begin_col[iat];
    int begin = this->iat2iwt_[iat];
    for(int i = 0; i < size; ++i)
    {
#ifdef __MPI
        indexes[i] = this->local2global_col(icol + i) - begin;
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
