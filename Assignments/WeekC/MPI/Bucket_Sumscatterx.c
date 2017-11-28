#include <mpi.h>

int Bucket_Sum_scatterx( double *buf, int *displs, double *recvbuf, MPI_Comm comm)
{
  int me, nprocs, left, right, inbucket, outbucket, sumscatter_tag=9999;
  MPI_Status status;
  MPI_Request request;

  MPI_Comm_rank( comm, &me );
  MPI_Comm_size( comm, &nprocs );
  
  left  = ( me - 1 + nprocs ) % nprocs;
  right = ( me + 1 )          % nprocs;

  inbucket  = right;
  outbucket = me;
  
  for ( step=1; step<nprocs; step++ ){
    // Post receive for incoming bucket
    MPI_Irecv( &recvbuf[ displs[ inbucket ] ], displs[ inbucket+1 ] - displs[ inbucket ],
	       MPI_DOUBLE, right, sumscatter_tag, comm, &request );

    // Send outgoing bucket
    MPI_Send( &buf[ dipls[ outbucket ] ],  displs[ inbucket+1 ] - displs[ inbucket ],
	      MPI_DOUBLE, left, sumscatter_tag, comm );

    // 
    MPI_Wait( &request, &status );

    for ( int i=displs[ inbucket ]; i< displs[ inbucket+1 ]; i++ )
      buf[ i ] += recvbuf[ i ];
    
    outbucket = inbucket;
    inbucket = ( inbucket+1 ) % nprocs      
  }
}
