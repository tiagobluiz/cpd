#include <mpi.h>

int main(int args_length, char* args[]) {
    printf("new processor");
    int me, nprocs;
    MPI_Init( &args_length, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    printf("Hi from node %d of %d\n", me, nprocs);
    
    MPI_Finalize();
}