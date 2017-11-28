#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include<immintrin.h>

#define dabs( x ) ( (x) < 0 ? -(x) : x )

double FLA_Clock();      // This is a routine for extracting elapsed
			 // time borrowed from the libflame library

/* MaxAbsDiff computes the maximum absolute difference over all
   corresponding elements of two matrices */
double MaxAbsDiff( int, int, double *, int, double *, int );

/* RandomMatrix overwrites a matrix with random values */
void RandomMatrix( int, int, double *, int );

void dgemm_( char *, char *,                 // transA, transB
	     int *, int *, int *,            // m, n, k
	     double *, double *, int *,      // alpha, A, ldA
	               double *, int *,      //        B, ldB
	     double *, double *, int * );    // beta,  C, ldC

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

/* GemmWRapper is a common interface to all the implementations we will 
   develop so we don't have to keep rewriting this driver routine. */
void GemmWrapper( int, int, int, double *, int, double *, int, double *, int );

int main(int argc, char *argv[])
{
  int
    m, n, k,
    ldA, ldB, ldC,
    size, first, last, inc,
    i, irep,
    nrepeats;

  double
    d_one = 1.0,
    dtime, dtime_best, 
    diff, maxdiff = 0.0, gflops;

  double
    *A, *B, *C, *Cold, *Cref;

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes m=n=k=first to last in increments
     of inc will be performed.  (Actually, we are going to go from
     largest to smallest since this seems to give more reliable 
     timings.  */
  printf( "%% enter first, last, inc:" );
  scanf( "%d%d%d", &first, &last, &inc );
  printf( "%% %d %d %d \n", first, last, inc );

  i = 1;
  for ( size=last; size>= first; size-=inc ){
    /* we will only time cases where all three matrices are square */
    m = n = k = size;
    ldA = ldB = ldC = size;

    /* Gflops performed */
    gflops = 2.0 * m * n * k * 1e-09;

    /* Allocate space for the matrices.  We will use five arrays:
       A will be the address where A is stored.   Addressed with alpha(i,j).
       B will be the address where B is stored.   Addressed with beta(i,j).
       C will be the address where C is stored.   Addressed with gamma(i,j).

       Now, we will compute C = A B + C with via routine GemmWrapper
       and also with a reference implementation.  Therefore, we will
       utilize two more arrays:
 
       Cold will be the address where the original matrix C is
       stored.  

       Cref will be the address where the result of computing C = A B
       + C computed with the reference implementation will be stored.

       Final note: we use a special routine _mm_malloc to allocate
       space that is aligned so that the starting address is a
       multiple of 64.  This makes it faster to load data into vector
       registers.
    */

    A = ( double * ) _mm_malloc( ldA * k * sizeof( double ), 64 );
    B = ( double * ) _mm_malloc( ldB * n * sizeof( double ), 64 );
    C = ( double * ) _mm_malloc( ldC * n * sizeof( double ), 64 );
    Cold = ( double * ) _mm_malloc( ldC * n * sizeof( double ), 64 );
    Cref = ( double * ) _mm_malloc( ldC * n * sizeof( double ), 64 );

    /* Generate random matrix A */
    RandomMatrix( m, k, A, ldA );

    /* Generate random matrix B */
    RandomMatrix( k, n, B, ldB );

    /* Generate random matrix Cold */
    RandomMatrix( m, n, Cold, ldC );
    
    /* Time reference implementation provided by the BLAS library
       routine dgemm (double precision general matrix-matrix
       multiplicationn */
    for ( irep=0; irep<nrepeats; irep++ ){
      
      /* Copy matrix Cold to Cref */
      memcpy( Cref, Cold, ldC * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute Cref = A B + Cref */
      dgemm_( "No transpose", "No transpose",
	      &m, &n, &k,
	      &d_one, A, &ldA,
	              B, &ldB,
	      &d_one, Cref, &ldC );

      /* stop clock */
      dtime = FLA_Clock() - dtime;

      /* record the best time so far */
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }
  
    printf( "data_ref( %d, 1:3 ) = [ %d %le %le];\n",
	    i, n, dtime_best, gflops/dtime_best );
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.

    /* Time GemmWrapper */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector Cold to C */
      memcpy( C, Cold, ldC * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute C = A B + C */
      GemmWrapper( m, n, k, A, ldA, B, ldB, C, ldC );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = MaxAbsDiff( m, n, C, ldC, Cref, ldC );
    maxdiff = ( diff > maxdiff ? diff : maxdiff );
    
    printf( "data_Gemm( %d, 1:4 ) = [ %d %le %le %le];\n",
	    i, n, dtime_best, diff, gflops/dtime_best  );
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.

    /* Free the buffers */
    _mm_free( A );
    _mm_free( B );
    _mm_free( C );
    _mm_free( Cold );
    _mm_free( Cref );

    i++;
  }

  printf( "\n \% Maximum difference between reference and your implementation: %le.\n", maxdiff );
  
  exit( 0 );
}
