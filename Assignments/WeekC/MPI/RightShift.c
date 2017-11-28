#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
  int nprocs, me, left, right, outgoing;
  MPI_Status status;
  
  MPI_Init( &argc, &argv );

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &me );
  left = ( me == 0 ? nprocs-1 : me-1 );
  right = ( me == nprocs-1 ? me+1 : 0 );

  outgoing = me;

  MPI_Send( &outgoing, 1, MPI_INT, right, 0, MPI_COMM_WORLD );
  MPI_Recv( &incoming, 1, MPI_INT, left,  0, MPI_COMM_WORLD, &status );
  
  printf( "me = %d, incoming = %d\n", me, incoming );

  MPI_Finalize();
}
