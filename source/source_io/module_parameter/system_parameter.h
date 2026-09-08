#ifndef SYSTEM_PARAMETER_H
#define SYSTEM_PARAMETER_H

#include <ctime>
#include <string>

// Derived process and workflow state. None of these values is user input.
struct System_para
{
    int myrank = 0;
    int nproc = 1;
    int nthread_per_proc = 1;
    std::time_t start_time = 0;

    int nlocal = 0;
    bool two_fermi = false;
    int npol = 1;
    bool domag = false;
    bool domag_z = false;
    bool search_pbc = true;

    std::string global_in_card = "INPUT";
    std::string global_in_stru = "STRU";
    std::string global_out_dir;
    std::string global_readin_dir;
    std::string log_file = "running_get_h0.log";
};

#endif
