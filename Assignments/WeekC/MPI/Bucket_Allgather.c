#include <mpi.h>

int Bucket_Allgatherx( int *buf, int *displs, MPI_Comm comm)
{
  int me, nprocs, left, right, inbucket, outbucket, allgather_tag=9999;
  MPI_Status status;
  MPI_Request request;

  MPI_Comm_rank( comm, &me );
  MPI_Comm_size( comm, &nprocs );
  
  left  = ( me - 1 + nprocs ) % nprocs;
  right = ( me + 1 )          % nprocs;

  inbucket  = left;
  outbucket = me;
  
  for ( int step=1; step<nprocs; step++ ){
    // Post receive for incoming bucket
    MPI_Irecv( &buf[ displs[ inbucket ] ], displs[ inbucket+1 ] - displs[ inbucket ],
	       MPI_INT, left, allgather_tag, comm, &request );

    // Send outgoing bucket
    MPI_Send( &buf[ displs[ outbucket ] ],  displs[ inbucket+1 ] - displs[ inbucket ],
	      MPI_INT, right, allgather_tag, comm );

    // 
    MPI_Wait( &request, &status );

    outbucket = inbucket;
    inbucket = ( inbucket-1+nprocs ) % nprocs;
  }
}
