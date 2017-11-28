#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
  int nprocs, me, dest, my_number=999, incoming=-999;
  MPI_Status status;

  MPI_Init( &argc, &argv );

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  if ( me == 0 ){
    printf( "Input an integer:" );
    scanf( "%d", &my_number );

    /* send to all other processes */
    for ( dest=1; dest<nprocs; dest++ )
      MPI_Send( &my_number, 1, MPI_INT, dest, 0, MPI_COMM_WORLD );
  }
  else
    MPI_Recv( &incoming, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
 

  printf( "me = %d, my_message = %d, incoming = %d\n", me, my_message, incoming );

  MPI_Finalize();
}
