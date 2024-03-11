#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* load the C headers of ira lib */
#include "iralib_interf.h"

/*
  This is a small example program in C, that shows how to declare and call routines
  from the shared library shlib_ira.so, that are defined in library_ira.f90

   WARNING: this code is NOT very robust to read different file types, empty lines, etc.
   USE AT OWN RISK

   Note the input file containing the structure is passed as argument, not stdin. Run as:

       c_program.x  input_file

*/


int main( int argc, char **argv ){

  char *fname;
  fname = argv[1];
  FILE * fp = fopen( fname, "r" );

  /* ======================== */
  /* read the first structure */
  /* ======================== */
  int nat1;
  int *typ1;
  double **coords1;

  /* read nat */
  fscanf(fp, "%d\n", &nat1 );

  /* allocate memory */
  typ1=malloc( sizeof(int)*nat1);
  double *data1;
  data1=malloc(sizeof(double)*3*nat1);
  coords1=malloc( sizeof(double)*nat1);
  int n = 0;
  for( int i=0; i< nat1; i++){
    coords1[i] = &data1[n];
    n+=3;
  }

  /* read line, scratch */
  fscanf( fp, "%*s\n" );

  /* read atoms */
  for( int i=0; i< nat1; i++ ){
    fscanf( fp, "%d %lf %lf %lf\n", &typ1[i], &coords1[i][0], &coords1[i][1], &coords1[i][2] );
  }


  /* ========================= */
  /* read the second structure */
  /* ========================= */
  int nat2;
  int *typ2;
  double **coords2;

  /* read nat */
  fscanf(fp, "%d\n", &nat2 );

  /* allocate memory */
  typ2=malloc( sizeof(int)*nat2);
  double *data2;
  data2=malloc(sizeof(double)*3*nat2);
  coords2=malloc( sizeof(double)*nat2);
  n = 0;
  for( int i=0; i< nat2; i++){
    coords2[i] = &data2[n];
    n+=3;
  }

  /* read line, scratch */
  fscanf( fp, "%*s\n" );

  /* read atoms */
  for( int i=0; i< nat2; i++ ){
    fscanf( fp, "%d %lf %lf %lf\n", &typ2[i], &coords2[i][0], &coords2[i][1], &coords2[i][2] );
  }


  printf( "%d\n", nat1);
  printf("original structure1\n");
  for ( int i=0; i< nat1; i++)
    {
      printf( "%d %lf %lf %lf \n", typ1[i], coords1[i][0], coords1[i][1], coords1[i][2]);
    }
  printf( "%d\n", nat2);
  printf("original structure2\n");
  for ( int i=0; i< nat2; i++)
    {
      printf( "%d %lf %lf %lf \n", typ2[i], coords2[i][0], coords2[i][1], coords2[i][2]);
    }



  /* form the candidate arrays */
  int *candidate1;
  int *candidate2;

  candidate1=malloc( sizeof(int)*nat1 );
  candidate2=malloc( sizeof(int)*nat2 );

  for( int i = 0; i<nat1; i++ ){
    candidate1[i]=0;
  }
  for( int i = 0; i<nat2; i++ ){
    candidate2[i]=0;
  }

  /* if structures equal number of atoms, do default candidates (send -1) */
  if( nat1 == nat2){
    candidate1[0]=-1;
    candidate2[0]=-1;
  }
  /* if structures non-equal number of atoms, take index 1 from struc1 and all idx froms truc2 */
  /* NOTE: the indices should be starting at 1 (fortran style) */
  else{
    candidate1[0]=1;
    for( int i = 0; i<nat2; i++){
      candidate2[i]=i+1;
    }
  }


  /* allocate arrays for result */
  double *data_rotation;
  double *translation;
  int *perm;
  double hd;

  data_rotation = malloc( sizeof(double)*9);
  translation = malloc( sizeof(double)*3);
  perm = malloc( sizeof(int)*nat2);

  double kmax_factor;
  kmax_factor = 1.8;

  int err;

  /* =========================================== */
  /* main call to lib_match from library_ira.f90 */
  /* =========================================== */
  libira_match( nat1, typ1, &coords1[0][0], candidate1, \
                nat2, typ2, &coords2[0][0], candidate2, \
                kmax_factor, &data_rotation, &translation, &perm, &hd, &err );

  if( err != 0){
    printf(" err code: %d\n", err);
    return err;
  }

  printf(" \n");
  printf( "final Hausdorff distance: %lf\n", hd );
  printf(" \n");

  /* ============================================================================ */
  /* the output arrays are all with 1d shape, you need to reshape them as needed! */
  /* ============================================================================ */


  printf("rotation matrix:\n");
  int m=0;
  for( int i=0; i<3; i++){
    for( int j=0;j<3;j++){
      printf( " %f ", data_rotation[m] );
      m+=1;
    }
    printf("\n");
  }

  printf(" \n");
  printf("translation vector:\n");
  printf( "%f %f %f\n", translation[0], translation[1], translation[2] );



  /* free the memory */
  free(candidate1);
  free(candidate2);
  free(typ1);
  free(typ2);
  free(coords1[0]);
  free(coords2[0]);
  free(coords1);
  free(coords2);
  free( data_rotation );
  free( translation );
  free( perm );

}
