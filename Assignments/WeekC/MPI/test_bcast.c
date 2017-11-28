#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void MST_Bcast( void *, int, MPI_Datatype, int, MPI_Comm );

#define root 2
#define local_count 4

int main(int argc, char** argv)
{
  int nprocs, me, token, *buf;
  MPI_Status status;
  
  /* Extract MPI runtime parameters and update argc and argv
     leaving the other runtime parameters. */
  MPI_Init( &argc, &argv );

  /* Extract the number of processes */
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  /* Extract the rank of this process */
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  /* allocate the buffer */
  buf = (int *) malloc( local_count * nprocs * sizeof( int ) );

  if ( me == root ) /* Initialize the buffer on the root */
    for ( int i=0; i<local_count*nprocs; i++ )
      buf[ i ] = i;
  else
    for ( int i=0; i<local_count*nprocs; i++ )
      buf[ i ] = me;

  MST_Bcast( buf, local_count*nprocs, MPI_INT, root, MPI_COMM_WORLD );

  /* Print out the buffer, one process at a time */
  if ( me > 0 )
    MPI_Recv( &token, 1, MPI_INT, me-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
    
  printf( "me = %d\n", me );
  for ( int i=0; i<local_count*nprocs; i++ )
    printf( "%d ", buf[ i ] );
  printf( "\n\n" );
  fflush( stdout );
  
  if ( me < nprocs-1 )
    MPI_Send( &token, 1, MPI_INT, me+1, 0, MPI_COMM_WORLD );

  free( buf );
  
  /* Clean up the MPI environment. */
  MPI_Finalize();
}
