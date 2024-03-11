

/* declare functions from library_ira.f90 */
void lib_cshda( int nat1, int *typ1, double *coords1, \
                int nat2, int *typ2, double *coords2, \
                double thr, int **found, double **dists);

void lib_cshda_pbc( int nat1, int *typ1, double *coords1, \
                    int nat2, int *typ2, double *coords2, double * lat2, \
                    double thr, int **found, double **dists);

void lib_match(int nat1, int *typ1, double *coords1, int *cand1,\
               int nat2, int *typ2, double *coords2, int *cand2, \
               double km_factor, double **rmat, double **tr, int **perm, double *hd, int *cerr);

void lib_get_version( char *string, int *date );


/* declare functions from library_sofi.f90 */
void lib_compute_all( int nat, int *typ, double *coords, double sym_thr, \
                      int *nmat, double **mat_data, int **perm_data, \
                      char **op_data, int **n_data, int **p_data,       \
                      double **ax_data, double **angle_data, double **dmax_data, char **pg, \
                      double **prin_ax, int *cerr );

void lib_get_symm_ops( int nat, int *typ, double *coords, double sym_thr, \
                       int *nmat, double **mat_data, int *cerr );

void lib_get_pg( int nmat, double *mat_data, char **pg, double **prin_ax, int verb, int *cerr);

void lib_unique_ax_angle( int nmat, double *mat_data, \
                          char **op_data, double **ax_data, double **angle_data, int *cerr);

void lib_analmat( double *mat, char **op, int *n, int *p, double **ax, double *angle, int *cerr);

void lib_get_perm( int nat, int *typ, double *coords, \
                   int nmat, double *mat_data, int **perm_data, double **dmax_data);

void lib_get_combos( int nat, int *typ, double *coords, int nmat, double *mat_data, \
                     int *nmat_out, double **mat_out);

void lib_try_mat( int nat, int *typ, double *coords, double *rmat, double *dh, int **perm);

void lib_construct_operation( char *op, double *axis, double angle, double **rmat, int *cerr);

void lib_mat_combos( int nmat, double *mat_data, int *nmat_out, double **mat_out);

void lib_ext_bfield( int nmat, double *mat_data, double *b_vec, \
                     int *nmat_out, double **mat_out);

void lib_get_err_msg( int ierr, char** msg );

void lib_matrix_distance( const double *mat1, const double *mat2, double *dist);
