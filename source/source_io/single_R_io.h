#ifndef SINGLE_R_IO_H
#define SINGLE_R_IO_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include <map>

namespace ModuleIO
{
    template <typename T>
    void output_single_R(std::ofstream& ofs,
        const std::map<size_t, std::map<size_t, T>>& XR,
        const double& sparse_threshold,
        const bool& binary,
        const Parallel_Orbitals& pv,
        const bool& reduce = true);
    
    // multi-processor reduce and output, only rank 0 write to file
    // not recommended for large computation systems
    template <typename T>
    void output_single_R_reduce(std::ofstream& ofs,
        const std::map<size_t, std::map<size_t, T>>& XR,
        const double& sparse_threshold,
        const bool& binary,
        const Parallel_Orbitals& pv);

    // single processor output without reduce, for system under 30000 atoms
    // OR multi-processor output without reduce, each processor outputs its own part
    // recommended for large computation systems
    // all data will be put together with other program
    template <typename T>
    void output_single_R_non_reduce_txt(std::ofstream& ofs,
        const std::map<size_t, std::map<size_t, T>>& XR,
        const double& sparse_threshold,
        const Parallel_Orbitals& pv);
    template <typename T>
    void output_single_R_non_reduce_binary(std::ofstream& ofs,
        const std::map<size_t, std::map<size_t, T>>& XR,
        const double& sparse_threshold,
        const Parallel_Orbitals& pv);
    
}

#endif
