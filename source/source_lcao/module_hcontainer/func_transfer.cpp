#include "hcontainer_funcs.h"
#include "transfer.h"

#ifdef __MPI
#include <mpi.h>

#include <vector>

namespace hamilt
{
namespace
{

template <typename TR>
void gather_values(const HContainer<TR>& hR_parallel,
                   HContainer<TR>* hR_serial,
                   int serial_rank)
{
    int my_rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    HTransSerial<TR>* serial_transfer = nullptr;
    if (my_rank == serial_rank)
    {
        serial_transfer = new HTransSerial<TR>(size, hR_serial);
    }
    HTransPara<TR> parallel_transfer(size, const_cast<HContainer<TR>*>(&hR_parallel));

    if (my_rank == serial_rank)
    {
        std::vector<int> indexes;
        serial_transfer->cal_ap_indexes(serial_rank, &indexes);
        parallel_transfer.receive_ap_indexes(serial_rank, indexes.data(), indexes.size());
        parallel_transfer.cal_orb_indexes(serial_rank, &indexes);
        serial_transfer->receive_orb_indexes(serial_rank, indexes.data(), indexes.size());

        for (int rank = 0; rank < size; ++rank)
        {
            if (rank != serial_rank)
            {
                serial_transfer->send_ap_indexes(rank);
            }
        }
        for (int rank = 0; rank < size; ++rank)
        {
            if (rank != serial_rank)
            {
                serial_transfer->receive_orb_indexes(rank);
            }
        }
    }
    else
    {
        parallel_transfer.receive_ap_indexes(serial_rank);
        parallel_transfer.send_orb_indexes(serial_rank);
    }

    long max_size = 0;
    std::vector<TR> received_values;
    if (my_rank == serial_rank)
    {
        max_size = serial_transfer->get_max_size();
        received_values.resize(max_size * size);
    }
    MPI_Bcast(&max_size, 1, MPI_LONG, serial_rank, MPI_COMM_WORLD);

    std::vector<TR> sent_values(max_size);
    parallel_transfer.pack_data(serial_rank, sent_values.data());
    MPI_Gather(sent_values.data(),
               max_size,
               MPITraits<TR>::datatype(),
               received_values.data(),
               max_size,
               MPITraits<TR>::datatype(),
               serial_rank,
               MPI_COMM_WORLD);

    if (my_rank == serial_rank)
    {
        for (int rank = 0; rank < size; ++rank)
        {
            serial_transfer->receive_data(rank, received_values.data() + rank * max_size);
        }
        delete serial_transfer;
    }
}

} // namespace

template <typename TR>
void transferSerials2Parallels(const HContainer<TR>& hR_serial,
                              HContainer<TR>* hR_parallel)
{
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    HTransSerial<TR> serial_transfer(size, const_cast<HContainer<TR>*>(&hR_serial));
    HTransPara<TR> parallel_transfer(size, hR_parallel);

    std::vector<int> sent_indexes;
    std::vector<int> received_indexes;
    std::vector<int> sendcounts(size), recvcounts(size), sdispls(size), rdispls(size);
    for (int rank = 0; rank < size; ++rank)
    {
        std::vector<int> indexes;
        serial_transfer.cal_ap_indexes(rank, &indexes);
        sendcounts[rank] = static_cast<int>(indexes.size());
        sdispls[rank] = static_cast<int>(sent_indexes.size());
        sent_indexes.insert(sent_indexes.end(), indexes.begin(), indexes.end());
    }
    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int received_count = 0;
    for (int rank = 0; rank < size; ++rank)
    {
        rdispls[rank] = received_count;
        received_count += recvcounts[rank];
    }
    received_indexes.resize(received_count);
    MPI_Alltoallv(sent_indexes.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPI_INT,
                  received_indexes.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPI_INT,
                  MPI_COMM_WORLD);

    sent_indexes.clear();
    for (int rank = 0; rank < size; ++rank)
    {
        parallel_transfer.receive_ap_indexes(rank,
                                             received_indexes.data() + rdispls[rank],
                                             recvcounts[rank]);
        std::vector<int> indexes;
        parallel_transfer.cal_orb_indexes(rank, &indexes);
        sendcounts[rank] = static_cast<int>(indexes.size());
        sdispls[rank] = static_cast<int>(sent_indexes.size());
        sent_indexes.insert(sent_indexes.end(), indexes.begin(), indexes.end());
    }
    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    received_count = 0;
    for (int rank = 0; rank < size; ++rank)
    {
        rdispls[rank] = received_count;
        received_count += recvcounts[rank];
    }
    received_indexes.resize(received_count);
    MPI_Alltoallv(sent_indexes.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPI_INT,
                  received_indexes.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPI_INT,
                  MPI_COMM_WORLD);
    for (int rank = 0; rank < size; ++rank)
    {
        serial_transfer.receive_orb_indexes(rank,
                                            received_indexes.data() + rdispls[rank],
                                            recvcounts[rank]);
    }

    serial_transfer.get_value_size(sendcounts.data());
    int sent_count = 0;
    for (int rank = 0; rank < size; ++rank)
    {
        sdispls[rank] = sent_count;
        sent_count += sendcounts[rank];
    }
    parallel_transfer.get_value_size(recvcounts.data());
    received_count = 0;
    for (int rank = 0; rank < size; ++rank)
    {
        rdispls[rank] = received_count;
        received_count += recvcounts[rank];
    }

    std::vector<TR> sent_values(sent_count);
    std::vector<TR> received_values(received_count);
    for (int rank = 0; rank < size; ++rank)
    {
        if (sendcounts[rank] > 0)
        {
            serial_transfer.pack_data(rank, sent_values.data() + sdispls[rank]);
        }
    }
    MPI_Alltoallv(sent_values.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPITraits<TR>::datatype(),
                  received_values.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPITraits<TR>::datatype(),
                  MPI_COMM_WORLD);
    for (int rank = 0; rank < size; ++rank)
    {
        if (recvcounts[rank] > 0)
        {
            parallel_transfer.receive_data(rank, received_values.data() + rdispls[rank]);
        }
    }
}

template <typename TR>
void gatherParallels(const HContainer<TR>& hR_parallel,
                     HContainer<TR>* hR_serial,
                     int serial_rank)
{
    int my_rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::vector<int> local_ijrs = hR_parallel.get_ijr_info();
    if (my_rank == serial_rank)
    {
        hR_serial->insert_ijrs(&local_ijrs);
        for (int rank = 0; rank < size; ++rank)
        {
            if (rank == serial_rank)
            {
                continue;
            }
            int count = 0;
            MPI_Recv(&count, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<int> remote_ijrs(count);
            MPI_Recv(remote_ijrs.data(), count, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            hR_serial->insert_ijrs(&remote_ijrs);
        }
        hR_serial->allocate();
    }
    else
    {
        const int count = static_cast<int>(local_ijrs.size());
        MPI_Send(&count, 1, MPI_INT, serial_rank, 0, MPI_COMM_WORLD);
        MPI_Send(local_ijrs.data(), count, MPI_INT, serial_rank, 0, MPI_COMM_WORLD);
    }

    gather_values(hR_parallel, hR_serial, serial_rank);
}

template void transferSerials2Parallels(const HContainer<double>&, HContainer<double>*);
template void transferSerials2Parallels(const HContainer<std::complex<double>>&,
                                        HContainer<std::complex<double>>*);
template void gatherParallels(const HContainer<double>&, HContainer<double>*, int);
template void gatherParallels(const HContainer<std::complex<double>>&,
                              HContainer<std::complex<double>>*,
                              int);

} // namespace hamilt
#endif
