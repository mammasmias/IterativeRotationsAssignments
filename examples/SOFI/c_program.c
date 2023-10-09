#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* load the C headers of ira lib */
#include "iralib_interf.h"

/* a small function to slice strings */
void slice(const char* str, char* result, size_t start, size_t end)
{
  strncpy(result, str + start, end - start);
}


int main( void ){

  int nat;
  int *typ;
  double *coords_data;
  double **coords;
  double sym_thr;


  /* ---------------------------------------------- */
  /* read an xyz file with sym_thr in second line */
  scanf( "%d\n", &nat);

  /* alolcate typ */
  typ=malloc( sizeof(int)*nat);

  /* allocate data space which is contigus in memory, after cut it into (3,nat) for coords */
  coords_data = malloc( sizeof(double)*3*nat);
  /* allocate coords */
  coords=malloc( sizeof(double)*nat );
  /* coords[i] points to a stride in data for each atom */
  int n=0;
  for( int i=0; i< nat; i++)
    {
      coords[i] = &coords_data[n];
      n+=3;
    }

  /* read symmetry threshold from second line */
  scanf( "%lf\n", &sym_thr);


  for( int i=0; i<nat; i++)
    {
      scanf("%d %lf %lf %lf\n", &typ[i], &coords[i][0], &coords[i][1], &coords[i][2] );
    }

  printf( "%d\n", nat);
  printf("\n");
  for ( int i=0; i< nat; i++)
    {
      printf( "%d %lf %lf %lf \n", typ[i], coords[i][0], coords[i][1], coords[i][2]);
    }
  /* ---------------------------------------------- */


  /* shift coords to chosen origin, in this case the geometrical center */
  double gc[3];
  gc[0] = 0.0; gc[1]=0.0; gc[2]=0.0;
  for (int i=0; i < nat; i++)
    {
      gc[0] = gc[0] + coords[i][0];
      gc[1] = gc[1] + coords[i][1];
      gc[2] = gc[2] + coords[i][2];
    }
  gc[0] = gc[0]/nat; gc[1] = gc[1]/nat; gc[2] = gc[2]/nat;

  /* shift coords */
  for (int i=0; i < nat; i++)
    {
      coords[i][0] = coords[i][0] - gc[0];
      coords[i][1] = coords[i][1] - gc[1];
      coords[i][2] = coords[i][2] - gc[2];
    }



  /* nmax value comes from sofi_tools.f90, is used just for initial allocation */
  int nmax=200;

  int nmat;
  double *mat_data;
  int *perm_data;
  char *op_data;
  int *n_data;
  int *p_data;
  double *ax_data;
  double *angle_data;
  double *dmax_data;
  char *pg;

  /* allocate space up to nmax for all arrays, all arrays are 1d and contiguous */
  /* if you want to reshape the output data into proper dim-arrays, it's up to you. */
  mat_data = malloc( sizeof(double)*9*nmax );
  perm_data = malloc( sizeof(int)*nat*nmax);
  n_data = malloc( sizeof(int)*nmax);
  p_data = malloc( sizeof(int)*nmax);
  ax_data = malloc( sizeof(double)*3*nmax);
  angle_data = malloc( sizeof(double)*nmax);
  dmax_data = malloc( sizeof(double)*nmax);
  op_data = malloc(sizeof(char)*2*nmax+1);
  pg = malloc(sizeof(char)*11);



  /* call SOFI compute_all */
  lib_compute_all( nat, typ, &coords[0][0], sym_thr,     \
                   &nmat, &mat_data, &perm_data,  \
                   &op_data, &n_data, &p_data,    \
                   &ax_data, &angle_data, &dmax_data, &pg );

  printf( "%d\n",nmat);

  printf( "all data:\n");
  int m=0;
  int mm=0;
  int mp=0;
  char this_op[2]="";
  for ( int i=0; i<nmat; i++)
  {
    printf( " idx: %d\n", i );

    /* op_data is a single string, should be sliced into nmat substrings of size-2 each*/
    slice( op_data, this_op, i*2, i*2+2 );
    printf( " operation: %s %d ^ %d \n", this_op, n_data[i], p_data[i]);
    printf( " angle: %f\n",angle_data[i]);

    /* the axes data is concatenated for all symm */
    printf( " axis: %lf %lf %lf \n",ax_data[m],ax_data[m+1],ax_data[m+2]);
    m+=3;

    /* the matrix data is also concatenated into 1D */
    printf(" matrix:\n");
    for (int k=0;k<3;k++){
      for (int l=0;l<3;l++){
        printf( "  %f  ", mat_data[mm]);
        mm+=1;
      }
      printf("\n");
    }
    printf( " dmax  %f\n",dmax_data[i] );

    /* permutations are concatenated */
    printf( " permutation of atoms:\n");
    for( int k=0;k<nat;k++){
      printf( " %d ", perm_data[mp]);
      mp += 1;
    }
    printf("\n");
    printf("\n");

  }

  printf( "PG is: %s\n",pg);

  /* free all the allocated memory */
  free( mat_data );
  free(perm_data);
  free(op_data);
  free(n_data);
  free(p_data);
  free(ax_data);
  free(angle_data);
  free(dmax_data);
  free(pg);

  free(typ);
  free(coords[0]);
  free(coords);

  }
