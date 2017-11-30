#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void MST_Sum( int *, int, int, int *, MPI_Comm );

#define root 2
#define local_count 4

int main(int argc, char** argv)
{
  int nprocs, me, token, *buf, *tempbuf;
  MPI_Status status;
  
  /* Extract MPI runtime parameters and update argc and argv
     leaving the other runtime parameters. */
  MPI_Init( &argc, &argv );

  /* Extract the number of processes */
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  /* Extract the rank of this process */
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  /* allocate the buffers */
  buf     = (int *) malloc( local_count * nprocs * sizeof( int ) );
  tempbuf = (int *) malloc( local_count * nprocs * sizeof( int ) );

  for ( int i=0; i<local_count*nprocs; i++ )
    buf[ i ] = me*100 + i;

  MST_Sum( buf, local_count*nprocs, root, tempbuf, MPI_COMM_WORLD );

  /* Print out the buffer on the root */
  if ( me == root ) {  
    printf( "me = %d\n", me );
    for ( int i=0; i<local_count*nprocs; i++ )
      printf( "%d ", buf[ i ] );
    printf( "\n\n" );
    fflush( stdout );

    printf( "correct answer:\n" );
    for ( int i=0; i<local_count*nprocs; i++ )
      printf( "%d ", nprocs *( nprocs-1)/2 * 100 + nprocs*i );
    printf( "\n\n" );
    fflush( stdout );
  }

  free( buf );
  
  /* Clean up the MPI environment. */
  MPI_Finalize();
}
