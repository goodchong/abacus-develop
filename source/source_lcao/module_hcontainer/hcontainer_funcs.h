#pragma once

#include "source_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{

#ifdef __MPI
/** Redistribute rank-local integral containers into the 2D orbital layout. */
template <typename TR>
void transferSerials2Parallels(const HContainer<TR>& hR_serial,
                              HContainer<TR>* hR_parallel);

/** Gather and sum a distributed H(R) container onto one rank. */
template <typename TR>
void gatherParallels(const HContainer<TR>& hR_parallel,
                     HContainer<TR>* hR_serial,
                     int serial_rank);
#endif

} // namespace hamilt
