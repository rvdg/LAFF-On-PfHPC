#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"

#include "cblas.h"

#define UPLO FLA_LOWER_TRIANGULAR

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

void Symv_unb_var2( FLA_Obj, FLA_Obj, FLA_Obj );
void Symv_unb_var4( FLA_Obj, FLA_Obj, FLA_Obj );

int main(int argc, char *argv[])
{
  int n, nfirst, nlast, ninc, i, ii, jj, irep,
    nrepeats;

  double
    dtime, dtime_best, 
    diff;

  FLA_Obj
    Aobj, xobj, yobj, yold, yref;
  
  /* Initialize FLAME. */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  for ( n=nfirst; n<= nlast; n+=ninc ){

    /* Allocate space for the matrices and vectors */
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Aobj );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &xobj );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &yobj );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &yold );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, n, &yref );

    /* Generate random matrix A, and vectors x, and y */
    FLA_Random_matrix( Aobj );
    FLA_Random_matrix( xobj );
    FLA_Random_matrix( yold );



    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( yold, yref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      /* Compute yref = A x + yref where A is symmetric stored in the
	 lower triangular part of array A, by calling FLA_Symv.  The
	 result ends up in yref, which we will consider to be the
	 correct result. */
      FLA_Symv( UPLO, FLA_ONE, Aobj, xobj, FLA_ONE, yref );
      
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best );
    fflush( stdout );

    /* Time your unblocked Variant 2 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( yold, yobj );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      //    FLA_Symv( UPLO, FLA_ONE, Aobj, xobj, FLA_ONE, yobj );
      Symv_unb_var2( Aobj, xobj, yobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( yobj, yref );
  
    printf( "data_unb_var2( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );

    /* Time your unblocked Variant 4 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( yold, yobj );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      //    FLA_Symv( UPLO, FLA_ONE, Aobj, xobj, FLA_ONE, yobj );
      Symv_unb_var4( Aobj, xobj, yobj );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( yobj, yref );
  
    printf( "data_unb_var4( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    fflush( stdout );
    
    FLA_Obj_free( &Aobj );
    FLA_Obj_free( &xobj );
    FLA_Obj_free( &yobj );
    FLA_Obj_free( &yref );
    FLA_Obj_free( &yold );

    i++;
  }
  FLA_Finalize( );

  exit( 0 );
}
