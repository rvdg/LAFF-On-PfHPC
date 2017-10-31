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

void daxpy_( int *,                          // n
	     double *, double *, int *,      // alpha, x, incx
	               double *, int * );    //        y, incy
	      

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void Axpy( int n, double alpha, double *x, int incx, double *y, int incy );

int main(int argc, char *argv[])
{
  int
    n,
    incx, incy,
    size, first, last, inc,
    i, irep,
    nrepeats;

  double
    d_one = 1.0,
    dtime, dtime_best, 
    diff;

  double
    *xp, *yp, *yoldp, *yrefp;

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
    n = size;
    incx = incy = 1;

    /* Allocate space for the vectors.  We will use 4 vectors:
       xp will be the address where x is stored.
       yp will be the address where y is stored.  

       Now, we will compute y = alpha *x + y; with the routine in Axpy
       and also with a reference implementation.  Therefore, we will
       utilize two more vectors:
 
       yoldp will be the address where the original vector y is
       stored. 

       yrefp will be the address where the result of computing y = alpha * x 
       + y computed with the reference implementation will be stored.
    */
    
    xp = ( double * ) malloc( incx * n * sizeof( double ) );
    yp = ( double * ) malloc( incy * n * sizeof( double ) );
    yoldp = ( double * ) malloc( incy * n * sizeof( double ) );
    yrefp = ( double * ) malloc( incy * n * sizeof( double ) );

    /* Generate random vector x */
    RandomMatrix( n, 1, xp, n );

    /* Generate random vector yold */
    RandomMatrix( n, 1, yoldp, n );

    
    /* Time reference implementation provided by the BLAS library
       routine dgemm (double precision general matrix-matrix
       multiplicationn */
    for ( irep=0; irep<nrepeats; irep++ ){
      
      /* Copy vector yold to yref */
      memcpy( yrefp, yoldp, incy * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute yref = x + yref */
      daxpy_( &n, 
	      &d_one, xp, &incx,
	              yrefp, &incy);

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

    /* Time Axpy */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      memcpy( yp, yoldp, incy * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute C = A B + C */
      Axpy( n, d_one, xp, incx, yp, incy);

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = MaxAbsDiff( n, 1, yp, n, yrefp, n);
    
    printf( "data_Axpy( %d, 1:3 ) = [ %d %le %le];\n",
	    i, n, dtime_best, diff  );
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.

    /* Free the buffers */
    free( xp );
    free( yp );
    free( yoldp );
    free( yrefp );

    i++;
  }

  exit( 0 );
}
