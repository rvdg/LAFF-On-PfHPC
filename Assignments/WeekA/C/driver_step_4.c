#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
	     double *, double *, int *,      // alpha, A, lda
	               double *, int *,      //        B, ldb
	     double *, double *, int * );    // beta,  C, ldc

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

#define A( i,j ) *( ap + (j)*lda + (i) )          // map A( i,j )    to array ap    in column-major order
#define B( i,j ) *( bp + (j)*ldb + (i) )          // map B( i,j )    to array bp    in column-major order
#define C( i,j ) *( cp + (j)*ldc + (i) )          // map C( i,j )    to array cp    in column-major order
#define Cold( i,j ) *( coldp + (j)*ldc + (i) )    // map Cold( i,j ) to array coldp in column-major order
#define Cref( i,j ) *( crefp + (j)*ldc + (i) )    // map Cref( i,j ) to array crefp in column-major order

void GemmIJSapdot4x1Opt( int, int, int, double *, int, double *, int, double *, int );

int main(int argc, char *argv[])
{
  int
    m, n, k,
    lda, ldb, ldc,
    size, first, last, inc,
    i, irep,
    nrepeats;

  double
    d_one = 1.0,
    dtime, dtime_best, 
    diff;

  double
    *ap, *bp, *cp, *coldp, *crefp;

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes m=n=k=first to last in increments 
     of inc will be performed.  */
  printf( "%% enter first, last, inc:" );
  scanf( "%d%d%d", &first, &last, &inc );
  printf( "%% %d %d %d \n", first, last, inc );

  i = 1;
  for ( size=first; size<= last; size+=inc ){
    /* we will only time cases where all three matrices are square */
    m = n = k = size;
    lda = ldb = ldc = size;

    /* Allocate space for the matrices.  We will use five arrays:
       ap will be the address where A is stored.   Addressed with A(i,j).
       bp will be the address where B is stored.   Addressed with B(i,j).
       cp will be the address where C is stored.   Addressed with C(i,j).

       Now, we will compute C = A B + C with the routine in GemmIJSapdot4x1Opt
       and also with a reference implementation.  Therefore, we will
       utilize two more arrays:
 
       coldp will be the address where the original matrix C is
       stored.  Addressed with Cold(i,j).

       crefp will be the address where the result of computing C = A B
       + C computed with the reference implementation will be stored.
       Addressed with Cref(i,j).
    */
    
    ap = ( double * ) malloc( lda * k * sizeof( double ) );
    bp = ( double * ) malloc( ldb * n * sizeof( double ) );
    cp = ( double * ) malloc( ldc * n * sizeof( double ) );
    coldp = ( double * ) malloc( ldc * n * sizeof( double ) );
    crefp = ( double * ) malloc( ldc * n * sizeof( double ) );

    /* Generate random matrix A */
    RandomMatrix( m, k, ap, lda );

    /* Generate random matrix B */
    RandomMatrix( k, n, bp, ldb );

    /* Generate random matrix Cold */
    RandomMatrix( m, n, coldp, ldc );
    
    /* Time reference implementation provided by the BLAS library
       routine dgemm (double precision general matrix-matrix
       multiplicationn */
    for ( irep=0; irep<nrepeats; irep++ ){
      
      /* Copy matrix Cold to Cref */
      memcpy( crefp, coldp, ldc * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute Cref = A B + Cref */
      dgemm_( "No transpose", "No transpose",
	      &m, &n, &k,
	      &d_one, ap, &lda,
	              bp, &ldb,
	      &d_one, crefp, &ldc );

      /* stop clock */
      dtime = FLA_Clock() - dtime;

      /* record the best time so far */
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }
  
    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n",
	    i, n, dtime_best );
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.

    /* Time GemmIJSapdot4x1Opt */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector Cold to C */
      memcpy( cp, coldp, ldc * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute C = A B + C */
      GemmIJSapdot4x1Opt( m, n, k, ap, lda, bp, ldb, cp, ldc );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = MaxAbsDiff( m, n, cp, ldc, crefp, ldc );
    
    printf( "data_GemmIJSapdot4x1Opt( %d, 1:3 ) = [ %d %le %le];\n",
	    i, n, dtime_best, diff  );
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.

    /* Free the buffers */
    free( ap );
    free( bp );
    free( cp );
    free( coldp );
    free( crefp );

    i++;
  }

  exit( 0 );
}
