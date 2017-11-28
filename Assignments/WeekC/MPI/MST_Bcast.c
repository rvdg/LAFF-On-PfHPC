#include <mpi.h>

int MST_Bcast_aux( void *, int, MPI_Datatype, int, MPI_Comm, int, int );

int MST_Bcast( void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
{
  int nprocs;

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  MST_Bcast_aux( buf, count, datatype, root, comm, 0, nprocs-1 );
}

int MST_Bcast_aux( void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                  int left, int right )
{
  int me, mid, new_root, bcast_tag=9999;
  MPI_Status status;
  
  if ( left != right ){
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    mid = ( left + right ) / 2;
    new_root = ( root <= mid ? right : left );
    
    if ( me == root )
      MPI_Send( buf, count, datatype, new_root, bcast_tag, comm );

    if ( me == new_root )
      MPI_Recv( buf, count, datatype, root, bcast_tag, comm, &status );

    if ( me <= mid )  // me is in left component
      MST_Bcast_aux( buf, count, datatype,
                    ( root <= mid ? root : new_root ), comm, left, mid );
    else              // me is in right component
      MST_Bcast_aux( buf, count, datatype,
                    ( root <= mid ? new_root : root ), comm, mid+1, right );
  }
}
