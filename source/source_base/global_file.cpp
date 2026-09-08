#include "global_file.h"

#include "fs_compat.h"
#include "global_variable.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <cerrno>
#include <cstdlib>
#include <iostream>

namespace ModuleBase
{
namespace Global_File
{

void make_h0_output_dir(const int rank,
                        const std::string& output_dir,
                        const std::string& log_file)
{
    int ready = 1;
    if (rank == 0)
    {
        const int result = ModuleBase::make_directory(output_dir);
        ready = result == 0 || errno == EEXIST;
    }
#ifdef __MPI
    MPI_Bcast(&ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (!ready)
    {
        std::cerr << "Cannot create H0 output directory: " << output_dir << std::endl;
        std::exit(1);
    }

    if (rank == 0)
    {
        GlobalV::ofs_running.open((output_dir + log_file).c_str());
        GlobalV::ofs_warning.open((output_dir + "warning.log").c_str());
    }
}

void close_all_log(const int rank)
{
    if (rank == 0)
    {
        if (GlobalV::ofs_running)
        {
            GlobalV::ofs_running.close();
        }
        if (GlobalV::ofs_warning)
        {
            GlobalV::ofs_warning.close();
        }
    }
}

} // namespace Global_File
} // namespace ModuleBase
