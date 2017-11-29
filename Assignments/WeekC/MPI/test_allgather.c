#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void Bucket_Allgatherx( int *, int *, MPI_Comm );

#define local_count 4

int main(int argc, char** argv)
{
  int nprocs, me, token;
  int *buf;
  int *displs;
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

  /* allocate the displacement vector */
  displs = (int *) malloc( (nprocs+1) * sizeof( int ) );

  displs[ 0 ] = 0;
  for ( int i=1; i<=nprocs; i++ )
    displs[ i ] = displs[ i-1 ] + local_count;
  
  /* initialize the different parts of the vector */
  for ( int i=0; i<local_count*nprocs; i++ )
    buf[ i ] = 0;
  
  for ( int i=displs[ me ]; i<displs[ me+1 ]; i++ )
      buf[ i ] = i;

  Bucket_Allgatherx( buf, displs, MPI_COMM_WORLD );

  /* Print out the buffer, one process at a time */
  if ( me > 0 )
    MPI_Recv( &token, 1, MPI_INT, me-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
    
  printf( "me = %d\n", me );
  for ( int i=0; i<displs[ nprocs ]; i++ )
    printf( "%d ", buf[ i ] );
  printf( "\n\n" );
  fflush( stdout );
  
  if ( me < nprocs-1 )
    MPI_Send( &token, 1, MPI_INT, me+1, 0, MPI_COMM_WORLD );

  free( buf );
  
  /* Clean up the MPI environment. */
  MPI_Finalize();
}
