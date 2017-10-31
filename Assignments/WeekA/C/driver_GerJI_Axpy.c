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

void dger_ ( int *, int *,                   // m, n
	     double *, double *, int *,      // alpha, x, incx
	               double *, int *,      //        y, incy
	               double *, int *);     //        A, lda

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void GerJI_Axpy( int m, int n, double *x, int incx, double *y, int incy,
                               double *A, int ldA );

int main(int argc, char *argv[])
{
  int
    m, n,
    incx, incy, lda,
    size, first, last, inc,
    i, irep,
    nrepeats;

  double
    d_one = 1.0,
    dtime, dtime_best, 
    diff;

  double
    *xp, *yp, *Ap, *Aoldp, *Arefp;

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
    m = n = size;
    incx = incy = 1;
    lda = size;

    /* Allocate space for the vectors and matrices. We will use 
       2 vectors & 3 matrices:
       xp will be the address where x is stored.
       yp will be the address where y is stored.  
       Ap will be the address where A is stored.  

       Now, we will compute A = alpha * x * y^T + A; with the routine 
       in GerJI_Axpy and also with a reference implementation. 
       Therefore, we will utilize two more matrices:
 
       Aoldp will be the address where the original matrix A is
       stored. 

       Arefp will be the address where the result of computing 
       A = alpha * x * y^T + A computed with the reference 
       implementation will be stored.
    */
    
    xp = ( double * ) malloc( incx * m * sizeof( double ) );
    yp = ( double * ) malloc( incy * n * sizeof( double ) );
    Ap = ( double * ) malloc( lda * n * sizeof( double ) );
    Aoldp = ( double * ) malloc( lda * n * sizeof( double ) );
    Arefp = ( double * ) malloc( lda * n * sizeof( double ) );

    /* Generate random vector x */
    RandomMatrix( m, 1, xp, m );

    /* Generate random vector y */
    RandomMatrix( n, 1, yp, n );

    /* Generate random matrix Aold */
    RandomMatrix( m, n, Aoldp, lda );

    
    /* Time reference implementation provided by the BLAS library
       routine dgemm (double precision general matrix-matrix
       multiplicationn */
    for ( irep=0; irep<nrepeats; irep++ ){
      
      /* Copy matrix Aold to Aref */
      memcpy( Arefp, Aoldp, lda * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute A = x * y^T + A */
      dger_( &m, &n, 
	     &d_one, xp,    &incx,
	             yp,    &incy,
                     Arefp, &lda);

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
      /* Copy matrix Aold to A  */
      memcpy( Ap, Aoldp, lda * n * sizeof( double ) );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute A = x * y^T + A */
      GerJI_Axpy( m, n, xp, incx, yp, incy, Ap, lda);

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    /*for(int row = 0; row < m; row++){
      for(int col = 0; col < n; col++){
        printf("%f\t", Ap[ (col)*lda + row ] );
      }
      printf("\n");
    }

    printf("Aref\n");

    for(int row = 0; row < m; row++){
      for(int col = 0; col < n; col++){
        printf("%f\t", Arefp[ (col)*lda + row ] );
      }
      printf("\n");
    }

    printf("Aold\n");

    for(int row = 0; row < m; row++){
      for(int col = 0; col < n; col++){
        printf("%f\t", Aoldp[ (col)*lda + row ] );
      }
      printf("\n");
    }

   printf("x\n");

    for(int row = 0; row < m; row++){      
        printf("%f\t", xp[ row*incx ] );
      
    }
    printf("\n");

printf("y\n");

    for(int row = 0; row < n; row++){      
        printf("%f\t", yp[ row*incy ] );
      
    }
    printf("\n");	
*/

    diff = MaxAbsDiff( m, n, Ap, lda, Arefp, lda);
    
    printf( "data_GerJI_Axpy( %d, 1:3 ) = [ %d %le %le];\n",
	    i, n, dtime_best, diff  );
    fflush( stdout );  // We flush the output buffer because otherwise
		       // it may throw the timings of a next
		       // experiment.

    /* Free the buffers */
    free( xp );
    free( yp );
    free( Ap );
    free( Aoldp );
    free( Arefp );

    i++;
  }

  exit( 0 );
}
