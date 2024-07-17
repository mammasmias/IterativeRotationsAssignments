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


  /* check IRA version */

  /* char *version; */
  /* int date; */
  /* version = malloc( sizeof(char)*10); */
  /* libira_get_version( &version[0], &date ); */
  /* printf( "version %s\n", version ); */
  /* printf( "date %d\n", date); */
  /* free( version ); */



  int nat;               //! number of atoms
  int *typ;              //!  integer of atomic types
  double *coords_data;   //! atomic positions
  double **coords;
  double sym_thr;        //! threshold for symmetry, sym_thr is read from second line of the input file

  // set symmetry threshold
  sym_thr=0.05;


  /* ---------------------------------------------- */
  /* read an xyz file with empty second line */
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

  /* read empty line */
  fscanf( stdin, "\n");

  /* read coords */
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
  double gc[3]; //! origin point (geometric center)

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


  int nmax=400;          //! max space for initial allocations, value from sofi_tools.f90
  int nmat;              //! total number of symmetry operations
  double *mat_data;      //! list of symmetry operations
  int *perm_data;        //! list of permutations
  char *op_data;         //! list of Op
  int *n_data;           //! list of n values
  int *p_data;           //! list of p values
  double *ax_data;       //! list of axes
  double *angle_data;    //! list of angles
  double *dHausdorff_data;     //! list of dHausdorff values
  char *pg;              //! point group
  int n_prin_ax;         //! number of principal axes
  double *prin_ax;       //! list of principal axes
  int cerr;
  int prescreen_ih;

  /* Allocate space up to `nmax` for all arrays, all arrays are 1d and contiguous. */
  /* On output from SOFI, only values up to `nmat` are filled. */
  /* If you want to reshape the output data into proper dim-arrays, it's up to you. */
  mat_data = malloc( sizeof(double)*9*nmax );
  perm_data = malloc( sizeof(int)*nat*nmax);
  n_data = malloc( sizeof(int)*nmax);
  p_data = malloc( sizeof(int)*nmax);
  ax_data = malloc( sizeof(double)*3*nmax);
  angle_data = malloc( sizeof(double)*nmax);
  dHausdorff_data = malloc( sizeof(double)*nmax);
  op_data = malloc(sizeof(char)*1*nmax+1);
  pg = malloc(sizeof(char)*11);
  prin_ax = malloc(sizeof(double)*3*nmax);


  prescreen_ih = 0;

  /* call SOFI compute_all */
  libira_compute_all( nat, typ, &coords[0][0], sym_thr, prescreen_ih, \
                      &nmat, &mat_data, &perm_data,  \
                      &op_data, &n_data, &p_data,    \
                      &ax_data, &angle_data, &dHausdorff_data, &pg, &n_prin_ax, &prin_ax, &cerr );

  if( cerr < 0 ){
    char* msg;
    msg = malloc( sizeof(char)*128);
    libira_get_err_msg( cerr, &msg );
    printf( "%s\n", msg);
    return cerr;
  }

  printf( "%d\n",nmat);

  printf( "all data:\n");
  int m=0;
  int mm=0;
  int mp=0;
  char this_op[1]="";
  for ( int i=0; i<nmat; i++)
  {
    printf( " idx: %d\n", i+1 );

    /* op_data is a single string, should be sliced into nmat substrings of size-1 each*/
    slice( op_data, this_op, i*1, i*1+1 );
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
    printf( " dHausdorff  %f\n",dHausdorff_data[i] );

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
  printf( "List of principal axes, N_prin_ax = %d :\n", n_prin_ax);
  m=0;
  for( int i=0; i < n_prin_ax; i++ ){
     printf( "Principal axis: %lf %lf %lf\n", prin_ax[m+0], prin_ax[m+1], prin_ax[m+2] );
     m+=3;
  }

  /* free all the allocated memory */
  free( mat_data );
  free(perm_data);
  free(op_data);
  free(n_data);
  free(p_data);
  free(ax_data);
  free(angle_data);
  free(dHausdorff_data);
  free(pg);
  free(prin_ax);

  free(typ);
  free(coords[0]);
  free(coords);


  }
