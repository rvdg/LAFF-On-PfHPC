#include <mpi.h>
int MST_Scatterx_aux( int *, int *, int, MPI_Comm, int, int );

int MST_Scatterx( int *buf, int *displs, int root, MPI_Comm comm)
{
  int nprocs;

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MST_Scatterx_aux( buf, displs, root, comm, 0, nprocs-1 );
}

int MST_Scatterx_aux( int *buf, int *displs, 
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
        MPI_Send( &buf[ displs[ mid+1 ] ], displs[ right+1 ] - displs[  mid+1 ],
                  MPI_INT, new_root, scat_tag, comm );
      else
        MPI_Send( &buf[ displs[ left ] ], displs[ mid+1 ] - displs[ left ],
                  MPI_INT, new_root, scat_tag, comm );
    
    if ( me == new_root )
      if ( me > mid )
	MPI_Recv( &buf[ displs[ mid+1 ] ], displs[ right+1 ] - displs[ mid+1 ],
                  MPI_INT, root, scat_tag, comm, &status );
      else
	MPI_Recv( &buf[ displs[ left ] ], displs[ mid+1 ] - displs[ left ],
                  MPI_INT, root, scat_tag, comm, &status );
    
    if ( me <= mid )  //me is in left component
      MST_Scatterx_aux( buf, displs,
                        ( root <= mid ? root : new_root ), comm, left, mid );
    else              // me is in right component
      MST_Scatterx_aux( buf, displs,
                        ( root <= mid ? new_root : root ), comm, mid+1, right );
  }
}
