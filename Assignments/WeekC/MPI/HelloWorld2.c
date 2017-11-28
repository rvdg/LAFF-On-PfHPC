#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
  int nprocs, me, left, right;

  /* Extract MPI runtime parameters and update argc and argv
     leaving the other runtime parameters. */
  MPI_Init( &argc, &argv );

  /* Extract the number of processes */
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  /* Extract the rank of this process */
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  /*  Hello! */
  printf( "Hello World from %d of %d processes\n", me, nprocs );

  /* Clean up the MPI environment. */
  MPI_Finalize();
}
