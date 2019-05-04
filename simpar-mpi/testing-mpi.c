#include "mpi.h"
#include <stdio.h>

struct Partstruct
{
    char c;
    double d;
    char b;
};

int main(int argc, char *argv[])
{
    struct Partstruct particle[1000];
    int i, j, myrank;
    MPI_Status status;
    MPI_Datatype Particletype;
    MPI_Datatype type[3] = { MPI_CHAR, MPI_DOUBLE, MPI_CHAR };
    int blocklen[3] = { 1, 6, 7 };
    MPI_Aint disp[3];

    MPI_Init(&argc, &argv);

    disp[0] = &particle[0].c - &particle[0];
    disp[1] = &particle[0].d - &particle[0];
    disp[2] = &particle[0].b - &particle[0];
    MPI_Type_create_struct(3, blocklen, disp, type, &Particletype);
    MPI_Type_commit(&Particletype);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        MPI_Send(particle, 1000, Particletype, 1, 123, MPI_COMM_WORLD);
    }
    else if (myrank == 1)
    {
        MPI_Recv(particle, 1000, Particletype, 0, 123, MPI_COMM_WORLD, &status);
    }
    MPI_Finalize();
    return 0;
}
