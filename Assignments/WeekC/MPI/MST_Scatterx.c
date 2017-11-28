#include <mpi.h>

int MST_Scatterx( char *buf, int *displs, int root, MPI_Comm comm)
{
  int nprocs;

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  MST_Scatterx_aux( buf, displs, root, comm, 0, nprocs-1 );
}

int MST_Scatterx_aux( char *buf, int *displs, 
		      int root, MPI_Comm comm, int left, int right )
{
  int me, mid, new_root, scat_tag=9999;
  MPI_Status status;
  
  if ( left != right ){
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    mid = ( left + right ) / 2;
    new_root = ( root <= mid ? right : left );
    
    if ( me == root )
      if ( me <= mid )
        MPI_Send( &buf[ mid+1 ], displs[ right+1 ] - displs[ mid+1],
                  MPI_CHAR, new_root, scat_tag, comm );
      else
        MPI_Send( &buf[ left ], displs[ mid+1 ] - displs[ left ],
                  MPI_CHAR, new_root, scat_tag, comm );
    
    if ( me == new_root )
      if ( me > mid )
        MPI_Recv( &buf[ mid+1 ], displs[ right+1 ] - displs[ mid+1],
                  MPI_CHAR, root, scat_tag, comm, &status );
      else
        MPI_Recv( &buf[ left ], displs[ mid+1 ] - displs[ left ],
                  MPI_CHAR, root, scat_tag, comm, &status );
	
    if ( me <= mid )  // me is in left component
      MST_Scatterx_aux( &buf[ left ], &displs[left], 
                        ( root <= mid ? root : new_root ), comm, left, mid );
    else              // me is in right component
      MST_Scatterx_aux( &buf[ mid+1 ], &displs[ mid+1], 
                        ( root <= mid ? new_root : root ), comm, mid+1, right );
  }
}
