!!
!! Copyright (C) 2023, MAMMASMIAS Consortium
!! Written by: Miha Gunde
!!
!! SPDX-License-Identifier: GPL-3.0-or-later
!! SPDX-License-Identifier: Apache-2.0
!! See the file LICENSE.txt for further information.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!

!!==============================================================================
!! This file defines the wrappers to SOFI routines, which are part of the
!! shared library libira.so.
!! The routines here are defined with "bind(C)" and "use iso_c_binding", and
!! they effectively behave as C-code. The arrays passed to these routines are
!! assumed to be already allocated by the caller, and are assumed to have
!! C-shape. Take care of transposed matrices, starting indices (1 or 0), etc.
!! Likewise, the output arrays are assumed to be already allocated by the caller,
!! and it's up to the caller to assure the correct shape and indices.
!! The routines here receive C-pointers to the memory where output should be
!! stored. Therefore, the output data appears as "intent(in)".
!!==============================================================================

!!

!> @details
!! Perform the full computation by SOFI.
!! This is a wrapper to sofi_compute_all() from sofi_routines.f90
!!
!! @note
!!  the ``nmax`` refers to the value in sofi_tools.f90, which ``nmax=200`` by default.
!!
!! The size of variables on input needs to be at least this much.
!! The actual output is written to the first ``n_mat`` elements of each corresponding array.
!!
!! @param[in] nat                 :: number of atoms
!! @param[in] typ(nat)            :: atomic types
!! @param[in] coords(3,nat)       :: atomic positions
!! @param[in] sym_thr             :: threshold for finding symmetries
!!                                   (not taken into account when making combinations)
!! @param[out] n_mat              :: number of symmetries found
!! @param[in] mat_list(3,3,nmax)  :: symmetry matrices in C order
!! @param[in] perm_list(nat,nmax) :: permutation of atoms after applying each symmetry, in C format (start at 0)
!! @param[in] op_list(nmax)       :: Character "Op" from the Schoenflies notation: Op n^p
!!                                   (E = identity, I = inversion, C = rotation, S = (roto-)reflection )
!! @param[in] n_list(nmax)        :: Schoenflies n value
!! @param[in] p_list(nmax)        :: Schoenflies p value
!! @param[in] ax_list(3,nmax)     :: axis of operation of each symmetry
!! @param[in] angle_list(nmax)    :: angle of each symmetry, in units of 1/2pi,
!!                                   i.e. angle=0.333 is 1/3 of full circle, or 120 degrees
!! @param[in] dmax_list(nmax)     :: max difference of atomic positions of before/after symm transformation
!! @param[in] pg                  :: name of Point group, e.g. D6h
!! @param[in] prin_ax(3)          :: principal axis of the PG
!! @param[out] cerr               :: error value, negative on error, zero otherwise
!!
!! @returns n_mat, mat_list, per_list, op_list, n_list, p_list, ax_list, angle_list, dmax_list, pg, prin_ax
!!
!! C-header:
!! ~~~~~~~~~~~~~~~{.c}
!! void libira_compute_all( int nat, int *typ, double *coords, double sym_thr, \
!!                       int *n_mat, double **mat_data, int **perm_data, \
!!                       char **op_data, int **n_data, int **p_data,       \
!!                       double **ax_data, double **angle_data, double **dmax_data, char **pg, \
!!                       double **prin_ax, int *cerr );
!! ~~~~~~~~~~~~~~~
!!
subroutine libira_compute_all( nat, typ, coords, sym_thr, &
                            n_mat, mat_list, perm_list, &
                            op_list, n_list, p_list, &
                            ax_list, angle_list, dmax_list, pg, prin_ax, &
                            cerr ) bind(C, name="libira_compute_all")
  use iso_c_binding
  use sofi_tools, only: nmax
  use err_module
  implicit none
  integer( c_int ), value, intent(in) :: nat
  type( c_ptr ), value, intent(in) :: typ
  type( c_ptr ), value, intent(in) :: coords
  real( c_double ), value, intent(in) :: sym_thr

  integer( c_int ), intent(out) :: n_mat
  type( c_ptr ), intent(in) :: mat_list
  type( c_ptr ), intent(in) :: perm_list
  type( c_ptr ), intent(in) :: op_list
  type( c_ptr ), intent(in) :: n_list
  type( c_ptr ), intent(in) :: p_list
  type( c_ptr ), intent(in) :: ax_list
  type( c_ptr ), intent(in) :: angle_list
  type( c_ptr ), intent(in) :: dmax_list
  type( c_ptr ), intent(in) :: pg
  type( c_ptr ), intent(in) :: prin_ax
  integer( c_int ), intent(out) :: cerr

  !! pointers for c input arrays
  integer(c_int), dimension(:), pointer :: ptyp
  real( c_double), dimension(:,:), pointer :: pcoords
  real( c_double ), dimension(:,:,:), pointer :: pmat_list
  integer( c_int ), dimension(:,:), pointer :: pperm_list
  integer( c_int ), dimension(:), pointer :: pn_list
  integer( c_int ), dimension(:), pointer :: pp_list
  real( c_double ), dimension(:,:), pointer :: pax_list
  real( c_double ), dimension(:), pointer :: pangle_list
  real( c_double ), dimension(:), pointer :: pdmax_list
  real( c_double ), dimension(:), pointer :: pprin_ax
  character(len=1, kind=c_char), dimension(:), pointer :: pg_char
  character(len=1, kind=c_char), dimension(:), pointer :: op_char

  !! some f-defined memory
  integer( c_int ) :: nb
  character(len=10) :: f_pg
  character(len=1), dimension(nmax) :: fop_list
  character(len=1) :: this_op

  integer :: i, m, n, lenc, ierr
  integer :: len_opstr

  ! write(*,*) "entering lib"
  !! local arrays for computation

  !!
  !! receive input
  !!
  call c_f_pointer( typ, ptyp, [nat] )
  call c_f_pointer( coords, pcoords, [3,nat] )
  !!
  !! set pointers from c to f
  !!
  call c_f_pointer( mat_list, pmat_list, [3,3,nmax] )
  call c_f_pointer( perm_list, pperm_list, [nat,nmax] )
  call c_f_pointer( n_list, pn_list, [nmax] )
  call c_f_pointer( p_list, pp_list, [nmax] )
  call c_f_pointer( ax_list, pax_list, [3,nmax] )
  call c_f_pointer( angle_list, pangle_list, [nmax] )
  call c_f_pointer( dmax_list, pdmax_list, [nmax] )
  call c_f_pointer( pg, pg_char, [10+1] )
  call c_f_pointer( prin_ax, pprin_ax, [3] )

  !!
  !! compute SOFI
  !!
  call sofi_compute_all( nat, ptyp, pcoords, sym_thr, &
       nb, pmat_list, pperm_list, &
       fop_list, pn_list, pp_list, &
       pax_list, pangle_list, pdmax_list, f_pg, pprin_ax, ierr )

  cerr = int( ierr, c_int )
  if( ierr /= 0 ) then
     !! i think here cant call a function like this... call to sofi_get_err_msg(..)
     write(*,*) get_err_msg( ierr )
     return
  end if


  !! set pg string
  lenc=len_trim(f_pg)
  do i=1, lenc
     pg_char(i) = f_pg(i:i)
  end do
  pg_char(lenc+1)=c_null_char

  !! compact the op_list into single string to pass to C
  len_opstr = len(this_op)
  lenc=nmax*len_opstr
  call c_f_pointer( op_list, op_char, [lenc+1] )
  n=1
  do i = 1, nb
     this_op=fop_list(i)
     do m = 1, len_opstr
        op_char(n) = this_op(m:m)
        n = n + 1
     end do
  end do
  op_char(n)=c_null_char

  !! set output value for nmat
  n_mat = nb

  !! transpose all output matrices to C-order
  do i = 1, nb
     pmat_list(:,:,i) = transpose( pmat_list(:,:,i))
  end do

  !! return C-style indices
  pperm_list = pperm_list-1

  ! write(*,*) "exiting lib"
end subroutine libira_compute_all


!> @details
!!
!! Get the list of symmetry operations.
!! This is a wrapper to sofi_get_symmops from sofi_routines.f90
!!
!! @param[in] nat                 :: number of atoms
!! @param[in] typ(nat)            :: atomic types
!! @param[in] coords(3,nat)       :: atomic positions
!! @param[in] sym_thr             :: threshold for finding symmetries
!!                                   (not taken into account when making combinations)
!! @param[in] n_mat               :: number of symmetries found
!! @param[in] mat_list(3,3,nmax)  :: symmetry matrices in C format
!! @returns n_mat, mat_list
!!
!! C header:
!!~~~~~~~~~~~~~~~{.c}
!! void libira_get_symm_ops( int nat, int *typ, double *coords, double symm_thr, \
!!                        int *n_mat, double **mat_list, int *cerr );
!!~~~~~~~~~~~~~~~
!!
subroutine libira_get_symm_ops(nat, typ, coords, symm_thr, n_mat, mat_list, cerr )&
     bind(C, name="libira_get_symm_ops")
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  !! "input" structure
  integer( c_int ), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  real( c_double ), value, intent(in) :: symm_thr
  !! "output"
  integer( c_int ), intent(out) :: n_mat
  type( c_ptr ),    intent(in) :: mat_list
  integer( c_int ), intent(out) :: cerr

  !! F pointers
  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  real( c_double), dimension(:,:,:), pointer :: ptr_op
  !! memory in f
  integer(c_int) :: n, i
  integer :: ierr

  ! write(*,*) 'enter symmop'
  n_mat = 0_c_int
  !!
  !! receive input
  !!
  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_coords, [3,nat] )

  call c_f_pointer( mat_list, ptr_op, [3,3,nmax] )

  call sofi_get_symmops( nat, ptr_typ, ptr_coords, symm_thr, n, ptr_op, ierr )
  cerr = int( ierr, c_int )
  if( ierr /= 0 ) then
     write(*,*) "at",__FILE__,"line:",__LINE__
     return
  end if


  !! set output data
  n_mat=n

  !! transpose output matrices to C-order
  do i = 1, n
     ptr_op(:,:,i) = transpose( ptr_op(:,:,i) )
  end do

end subroutine libira_get_symm_ops


!> @details
!!
!! Get the point group from a list of matrices, and its principal axis.
!! This is a wrapper to sofi_get_pg() from sofi_routines.f90
!!
!! @param[in] n_mat :: number of matrices
!! @param[in] cptr_op_list(3,3,n_mat) :: list of matrices in C order
!! @param[in] ppg :: point group
!! @param[in] px :: principal axis
!! @param[in] verbose :: flag of verbosity
!! @returns ppg, px
!!
!! C-header:
!!~~~~~~~~~~~~~~{.c}
!! void libira_get_pg( int n_mat, double *mat_data, char **pg, double **prin_ax, int verb, int *cerr);
!!~~~~~~~~~~~~~~
!!
subroutine libira_get_pg( n_mat, cptr_op_list, ppg, px, verbose, cerr )bind(C, name="libira_get_pg")
  use iso_c_binding
  implicit none
  integer( c_int ), value, intent(in) :: n_mat
  type( c_ptr ), value, intent(in) :: cptr_op_list
  !! c ptr to write output
  type( c_ptr ), intent(in) :: ppg
  type( c_ptr ), intent(in) :: px
  logical( c_bool ), value, intent(in) :: verbose
  integer( c_int ), intent(out) :: cerr

  !! f pointers
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real(c_double), dimension(:,:,:), pointer :: op_list
  ! character(len=10, kind=c_char), pointer :: str
  character(len=1, kind=c_char), dimension(:), pointer :: pg_char
  real( c_double ), dimension(:), pointer :: prin_ax
  !! local
  integer :: i, n, ierr
  character(len=10) :: pg
  logical :: verb

  !! receive input
  call c_f_pointer( cptr_op_list, p_lvl2, [9,n_mat])


  !! put the array into proper fortran shape
  allocate( op_list(1:3,1:3,1:n_mat))
  do n = 1, n_mat
     op_list(:,:,n)=reshape( p_lvl2(:,n),[3,3])
     !! transpose because it is C-style on input
     op_list(:,:,n) = transpose( op_list(:,:,n) )
  end do

  call c_f_pointer( ppg, pg_char, [11] )
  call c_f_pointer( px, prin_ax, [3] )

  verb = verbose
  call sofi_get_pg( n_mat, op_list, pg, prin_ax, verb, ierr )
  cerr = int( ierr, c_int )
  if( ierr /= 0 ) then
     write(*,*) "at",__FILE__,"line:",__LINE__
     return
  end if


  n = len_trim(pg)
  do i = 1, n
     pg_char(i)=pg(i:i)
  end do
  pg_char(n+1)=c_null_char

end subroutine libira_get_pg


subroutine libira_unique_ax_angle( n_mat, cptr_mat_list, op_out, ax_out, angle_out, cerr ) &
     bind(C,name="libira_unique_ax_angle")
  use iso_c_binding
  use err_module
  implicit none
  integer(c_int), value, intent(in) :: n_mat
  type( c_ptr ), value, intent(in) :: cptr_mat_list

  type( c_ptr ), intent(in) :: op_out
  type( c_ptr ), intent(in) :: ax_out
  type( c_ptr ), intent(in) :: angle_out
  integer( c_int ), intent(out) :: cerr

  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real(c_double), dimension(:,:,:), pointer :: mat_list
  real( c_double ), dimension(:,:), pointer :: ptr_ax
  real( c_double ), dimension(:), pointer :: ptr_angle
  character(len=1, kind=c_char), dimension(:), pointer :: op_char

  integer :: i, n, m, len_opstr, ierr
  character(len=1), dimension(n_mat) :: fop_list
  character(len=1) :: this_op

  !! receive input
  call c_f_pointer( cptr_mat_list, p_lvl2, [9,n_mat] )
  call c_f_pointer( ax_out, ptr_ax, [3,n_mat] )
  call c_f_pointer( angle_out, ptr_angle, [n_mat] )


  !! put arrays into proder fortran order
  allocate( mat_list(1:3,1:3,1:n_mat))
  do n = 1, n_mat
     mat_list(:,:,n)=reshape( p_lvl2(:,n),[3,3])
     !! transpose becuse it's C-order on input
     mat_list(:,:,n) = transpose(mat_list(:,:,n) )
  end do

  call sofi_unique_ax_angle( n_mat, mat_list, fop_list, ptr_ax, ptr_angle, ierr )
  cerr = int( ierr, c_int )
  if( ierr /= 0 ) then
     write(*,*) get_err_msg( ierr )
     return
  end if

  !! concatenate op strings into single string len(2*nmat+1)
  len_opstr = len(this_op)
  n=n_mat*len_opstr
  call c_f_pointer( op_out, op_char, [n+1] )
  n=1
  do i = 1, n_mat
     this_op = fop_list(i)
     do m = 1, len_opstr
        op_char(n) = this_op(m:m)
        n = n + 1
     end do
  end do
  op_char(n) = c_null_char

end subroutine libira_unique_ax_angle


!> @details
!!
!! Analyse the input 3x3 matrix, obtain Op n^p, axis, and angle.
!! This is a wrapper to sofi_analmat() from sofi_routines.f90
!!
!! @param[in] c_rmat(3,3) :: 3x3 input matrix in C order
!! @param[in] c_op :: Schoenflies symbol of operation E, I, C, S
!! @param[out] n :: order of operation
!! @param[out] p :: power
!! @param[in] c_ax(3) :: axis
!! @param[out] angle :: angle in units 1/(2pi), e.g. angle=0.5 is half circle
!! @param[out] cerr :: error value, zero on normal execution, negative otherwise
!! @returns c_op, n, p, c_ax, angle, cerr
!!
!! C-header:
!!~~~~~~~~~~~~~~~{.c}
!! void libira_analmat( double *mat, char **op, int *n, int *p, double **ax, double *angle, int *cerr);
!!~~~~~~~~~~~~~~~
!!
subroutine libira_analmat( c_rmat, c_op, n, p, c_ax, angle, cerr )bind(C,name="libira_analmat")
  use iso_c_binding
  implicit none
  !! in
  type( c_ptr ), value, intent(in) :: c_rmat
  !! out
  type( c_ptr ), intent(in) :: c_op
  integer(c_int), intent(out) :: n
  integer(c_int), intent(out) :: p
  type( c_ptr ), intent(in) :: c_ax
  real( c_double ), intent(out) :: angle
  integer( c_int ), intent(out) :: cerr

  !! f pointers
  real(c_double), dimension(:,:), pointer :: rmat
  real( c_double), dimension(:), pointer :: ax
  character(len=1, kind=c_char), dimension(:), pointer :: op_fptr

  !! local mem
  integer :: i, m, ierr
  character(len=1) :: op

  !! connect c to f
  call c_f_pointer( c_rmat, rmat, [3,3] )
  rmat = transpose(rmat)
  call c_f_pointer( c_op, op_fptr, [2] )
  call c_f_pointer( c_ax, ax, [3] )

  call sofi_analmat( rmat, op, n, p, ax, angle, ierr )
  cerr = int(ierr, c_int )
  if( ierr /= 0 )then
     write(*,*) "at:",__FILE__,"line:",__LINE__
     return
  end if

  m=len_trim(op)
  do i = 1, m
     op_fptr(i) = op(i:i)
  end do
  op_fptr(m+1)=c_null_char

end subroutine libira_analmat


subroutine libira_ext_bfield( n_mat, cop_list, cb_field, n_out, cop_out )&
     bind(C, name="libira_ext_bfield")
  use iso_c_binding
  implicit none

  !! input
  integer(c_int), value, intent(in) :: n_mat
  type( c_ptr ), value, intent(in) :: cop_list
  type( c_ptr ), value, intent(in) :: cb_field
  !!
  integer( c_int ), intent(out) :: n_out
  type( c_ptr ), intent(in) :: cop_out
  ! type( c_ptr ), value :: cop_out

  !! f ptrs
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double ), dimension(:), pointer :: b_field

  !! f memory
  real(c_double), dimension(:,:,:), pointer :: pmat_list
  real(c_double), dimension(3,3,n_mat) :: mat_list
  integer :: n
  integer :: m, i

  m=n_mat

  !! connect c to f
  call c_f_pointer( cop_list, p_lvl2, [9,n_mat] )

  call c_f_pointer( cb_field, b_field, [3] )
  call c_f_pointer( cop_out, pmat_list, [3,3,n_mat])

  ! allocate( mat_list(1:3,1:3,1:n_mat))
  do n = 1, n_mat
     mat_list(:,:,n)=reshape( p_lvl2(:,n),[3,3])
     mat_list(:,:,n) = transpose( mat_list(:,:,n) )
  end do

  call sofi_ext_Bfield( m, mat_list, b_field )

  !! re-transpose to C-style order
  do i = 1, m
     mat_list(:,:,i) = transpose( mat_list(:,:,i) )
  end do
  !! set data to ptr
  pmat_list = mat_list
  n_out = m



end subroutine libira_ext_bfield


!> @details Obtain the permutations and dmax for a list of matrices.
!! This is a wrapper to sofi_get_perm() from sofi_routines.f90
!!
!! mat_list in C order
!! perm_list in C order (index start at 0)
!!
!! C-header:
!!~~~~~~~~~~~~~~~~~~~~~~~{.c}
!! void libira_get_perm( int nat, int *typ, double *coords, \
!!                       int nmat, double *mat_data, int **perm_data, double **dmax_data);
!!~~~~~~~~~~~~~~~~~~~~~~~
!! @param[in] nat                 :: number of atoms
!! @param[in] typ(nat)            :: atomic types
!! @param[in] coords(3,nat)       :: atomic positions
!! @param[in] n_mat :: number of matrices on list
!! @param[in] mat_list(3,3,n_mat) :: list of matrices in C order
!! @param[in] perm_list(nat,n_mat) :: list of permutations for each matrix, in C order (start at 0)
!! @param[in] dmax_list(n_mat) :: list of dmax values for each matrix
subroutine libira_get_perm( nat, typ, coords, n_mat, mat_list, perm_list, dmax_list)&
     bind(C,name="libira_get_perm")
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  integer( c_int ), value, intent(in) :: n_mat
  type( c_ptr ), value,  intent(in) :: mat_list
  !! "output"
  type( c_ptr ), intent(in) :: perm_list
  type( c_ptr ), intent(in) :: dmax_list

  !! F pointers
  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double), dimension(:,:,:), pointer :: ptr_op
  integer( c_int ), dimension(:,:), pointer :: pperm_list
  real( c_double ), dimension(:), pointer :: pdmax_list

  integer :: i

  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_coords, [3,nat] )

  call c_f_pointer( perm_list, pperm_list, [nat,n_mat] )
  call c_f_pointer( dmax_list, pdmax_list, [n_mat] )

  call c_f_pointer( mat_list, p_lvl2, [9,n_mat] )

  !! receive input array and reshape into fortran order
  allocate( ptr_op(1:3,1:3,1:n_mat))
  do i = 1, n_mat
     ptr_op(:,:,i) = reshape( p_lvl2(:,i),[3,3])
     !! transpose because it is C-order on input
     ptr_op(:,:,i) = transpose( ptr_op(:,:,i) )
  end do

  call sofi_get_perm( nat, ptr_typ, ptr_coords, n_mat, ptr_op, pperm_list, pdmax_list )

  !! output C-style: start at 0 indices
  pperm_list = pperm_list - 1

end subroutine libira_get_perm

!> @details perform matrix combinations with checking each matrix against the atomic structure.
!! This is a wrapper to sofi_get_combos() from sofi_routines.f90
!!
!! C-header:
!!~~~~~~~~~~~~~~~~~~~~{.c}
!! void libira_get_combos( int nat, int *typ, double *coords, int nmat, double *mat_data, \
!!                         int *nmat_out, double **mat_out);
!!~~~~~~~~~~~~~~~~~~~~
!! @param[in] nat                 :: number of atoms
!! @param[in] typ(nat)            :: atomic types
!! @param[in] coords(3,nat)       :: atomic positions
!! @param[in] n_mat_in :: number of input matrices
!! @param[in] mat_data(3,3,n_mat_in) :: list of input matrices in C order
!! @param[out] n_mat_out :: number of output matrices
!! @param[in] mat_out(3,3,n_mat_out) :: list of output matrices in C order
subroutine libira_get_combos( nat, typ, coords, n_mat_in, mat_data, n_mat_out, mat_out )&
     bind(C,name="libira_get_combos")
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  integer(c_int), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  integer( c_int ), value, intent(in) :: n_mat_in
  type( c_ptr ), value,  intent(in) :: mat_data
  integer( c_int ), intent(out) :: n_mat_out
  type( c_ptr ), intent(in) :: mat_out


  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double), dimension(:,:,:), pointer :: pcombo_out

  real(c_double), dimension(3,3,nmax) :: mat_list
  integer :: i, m

  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_coords, [3,nat] )

  call c_f_pointer( mat_data, p_lvl2, [9,n_mat_in] )

  do i = 1, n_mat_in
     mat_list(:,:,i) = reshape( p_lvl2(:,i),[3,3])
     mat_list(:,:,i) = transpose(mat_list(:,:,i))
  end do

  m = n_mat_in

  call sofi_get_combos( nat, ptr_typ, ptr_coords, m, mat_list )

  n_mat_out = m
  call c_f_pointer( mat_out, pcombo_out, [3,3,m] )
  pcombo_out = mat_list(1:3,1:3,1:m)


end subroutine libira_get_combos



!> @details
!! test matrix rmat: manually rotate structure and then compute cshda.
!!
!! C header:
!!~~~~~~~~~~~~~~~~{.c}
!! void libira_try_mat( int nat, int *typ, double *coords, double *rmat, double *dh, int **perm);
!!~~~~~~~~~~~~~~~~
!! @param[in] nat             :: number of atoms
!! @param[in] typ(nat)        :: atomic types
!! @param[in] coords(3,nat)   :: atomic positions
!! @param[in] rmat(3,3)       :: matrix to test in C order
!! @param[out] dh :: output Hausdorff from CShDA
!! @param[out] perm :: permutation of atomic indices after application of rmat in C order (start at 0)
subroutine libira_try_mat( nat, typ, coords, rmat, dh, perm )bind(C,name="libira_try_mat")
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  type( c_ptr ), value,  intent(in) :: rmat
  real(c_double), intent(out) :: dh
  type( c_ptr ), intent(in) :: perm

  !! f pointers
  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  integer( c_int ), dimension(:), pointer :: p_found

  integer :: i
  real( c_double), dimension(3,nat) :: tf_coords
  integer(c_int), dimension(nat) :: tf_typ
  integer( c_int), dimension(nat) :: found
  real(c_double),dimension(nat) :: dists
  real(c_double), dimension(3,3) :: frmat

  !! connect pointers
  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_coords, [3,nat] )

  call c_f_pointer( rmat, p_lvl2, [3,3] )
  frmat = transpose( p_lvl2 )

  !! rotate with rmat
  do i = 1, nat
     tf_coords(:,i) = matmul(frmat, ptr_coords(:,i))
     tf_typ(i) = ptr_typ(i)
  end do

  !! compute dist with cshda
  call cshda( nat, ptr_typ, ptr_coords, nat, tf_typ, tf_coords, 10.0, found, dists )

  ! do i = 1, nat
  !    write(*,*) i, found(i), dists(i)
  ! end do

  ! write(*,*) nat
  ! write(*,*) 'orig'
  ! do i = 1, nat
  !    write(*,*) ptr_typ(i), ptr_coords(:,i)
  ! end do
  ! write(*,*) nat
  ! write(*,*) 'tf'
  ! do i = 1, nat
  !    write(*,*) tf_typ(i), tf_coords(:,found(i))
  ! end do

  call c_f_pointer( perm, p_found, [nat] )
  p_found(:) = found(:) - 1

  dh = maxval(dists)


end subroutine libira_try_mat


!> @details
!! Construct a 3x3 matrix from input arguments Op, axis, angle
!! This is a wrapper to sofi_construct_operation() from sofi_routines.f90
!!
!! @param[in] op :: Schoenflies symbol E, I, C, S
!! @param[in] axis(3) :: desired axis
!! @param[in] angle :: desired angle in units 1/(2pi), e.g. angle=0.5 means half circle
!! @param[in] matrix(3,3) :: output matrix in C order
!! @returns matrix
!!
!! C-header:
!!~~~~~~~~~~~~~{.c}
!! void libira_construct_operation( char *op, double *axis, double angle, double **rmat, int *cerr);
!!~~~~~~~~~~~~~
!!
subroutine libira_construct_operation( op, axis, angle, matrix, cerr )bind(C,name="libira_construct_operation")
  use iso_c_binding
  implicit none
  interface
     FUNCTION c_strlen(str) BIND(C, name='strlen')
       IMPORT :: c_ptr, c_size_t
       IMPLICIT NONE
       TYPE(c_ptr), INTENT(IN), VALUE :: str
       INTEGER(c_size_t) :: c_strlen
     END FUNCTION c_strlen
  end interface

  type( c_ptr ), value, intent(in) :: op
  type( c_ptr ), value, intent(in) :: axis
  real( c_double ), value, intent(in) :: angle
  type( c_ptr ), intent(in) :: matrix
  integer( c_int ), intent(out) :: cerr

  real( c_double ), dimension(:), pointer :: ptr_ax
  real(c_double), dimension(:,:), pointer :: ptr_matrix

  real( c_double ), dimension(3,3) :: f_matrix
  character(len=1,kind=c_char), dimension(:), pointer :: ptr_op
  character(len=1) :: f_op
  integer( c_int ) :: i, n
  integer :: ierr

  n=c_strlen(op)
  call c_f_pointer( op, ptr_op, [n] )
  if( n .gt. 1 ) then
     write(*,*) "expected len-1 string as op, got:",n, ptr_op
     return
  end if

  f_op=""
  do i = 1, n
     f_op(i:i)=ptr_op(i)
  end do

  f_op=trim(adjustl(f_op))

  call c_f_pointer( axis, ptr_ax, [3] )
  call c_f_pointer( matrix, ptr_matrix, [3,3] )

  call sofi_construct_operation( f_op, ptr_ax, angle, f_matrix, ierr )
  cerr = int( ierr, c_int )
  if( ierr /= 0 ) then
     write(*,*) "at",__FILE__, "line:",__LINE__
     return
  end if


  !! return C-order matrix
  ptr_matrix = transpose(f_matrix)

end subroutine libira_construct_operation


!> @details
!! Perform matrix combinations until group completeness, without atomic structure.
!! This is a wrapper to sofi_mat_combos() from sofi_routines.f90
!!
!! C-header:
!!~~~~~~~~~~~~~~~~~~~~~~{.c}
!! void libira_mat_combos( int nmat, double *mat_data, int *nmat_out, double **mat_out);
!!~~~~~~~~~~~~~~~~~~~~~~
!!
!! @param[in] n_mat_in :: number of input matrices
!! @param[in] mat_data(3,3,n_mat_in) :: list of input matrices in C order
!! @param[out] n_mat_out :: number of output matrices
!! @param[in] mat_out(3,3,n_mat_out) :: list of output matrices in C order
subroutine libira_mat_combos( n_mat_in, mat_data, n_mat_out, mat_out )bind(C,name="libira_mat_combos")
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  integer( c_int ), value, intent(in) :: n_mat_in
  type( c_ptr ), value,  intent(in) :: mat_data
  integer( c_int ), intent(out) :: n_mat_out
  type( c_ptr ), intent(in) :: mat_out


  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double), dimension(:,:,:), pointer :: pcombo_out

  real(c_double), dimension(3,3,nmax) :: mat_list
  real(c_double), dimension(3,3,n_mat_in) :: mat_datap
  integer :: i, m

  call c_f_pointer( mat_data, p_lvl2, [9,n_mat_in] )

  do i = 1, n_mat_in
     mat_list(:,:,i) = reshape( p_lvl2(:,i),[3,3])
     mat_list(:,:,i) = transpose(mat_list(:,:,i))
  end do

  m = n_mat_in

  mat_datap(:,:,:) = mat_list(:,:,1:n_mat_in)

  call sofi_mat_combos( n_mat_in, mat_datap, m, mat_list )

  n_mat_out = m
  call c_f_pointer( mat_out, pcombo_out, [3,3,m] )
  pcombo_out = mat_list(1:3,1:3,1:m)


end subroutine libira_mat_combos

!> @details
!! Compute the distance between two matrices, using the matrix_distance function used
!! internally by SOFI, to determine if matrices are equal or not.
!!
!! C header:
!!~~~~~~~~~~~~~~~~~~~~~~{.c}
!! void libira_matrix_distance( double *mat1, double *mat2, double *dist);
!!~~~~~~~~~~~~~~~~~~~~~~
!! @param[in] mat1(3,3) :: matrix in C order
!! @param[in] mat2(3,3) :: matrix in C order
!! @param[out] dist :: distance
subroutine libira_matrix_distance( mat1, mat2, dist )bind(C,name="libira_matrix_distance" )
  use iso_c_binding
  use sofi_tools, only: matrix_distance
  implicit none
  real( c_double ), dimension(3,3), intent(in) :: mat1
  real( c_double ), dimension(3,3), intent(in) :: mat2
  real( c_double ), intent(out) :: dist

  real, dimension(3,3) :: frmat1, frmat2
  real :: fdist

  !! transpose from C input
  frmat1 = real( mat1 )
  frmat2 = real( mat2 )
  frmat1 = transpose( frmat1 )
  frmat2 = transpose( frmat2 )

  call matrix_distance( frmat1, frmat2, fdist )
  dist = real( fdist, c_double )

end subroutine libira_matrix_distance


!> @details
!! get the error message from IRA/SOFI associated to cerr value.
!! This is a wrapper to sofi_get_err_msg() from sofi_routines.f90
!!
!! @param[in] cerr :: integer error value
!! @param[in] cmsg :: error message string
!! @returns cmsg
!!
!! C-header:
!!~~~~~~~~~~~~~{.c}
!! void libira_get_err_msg( int ierr, char** msg );
!!~~~~~~~~~~~~~
!!
subroutine libira_get_err_msg( cerr, cmsg )bind(C, name="libira_get_err_msg")
  use iso_c_binding
  implicit none
  integer( c_int ), value, intent(in) :: cerr
  type( c_ptr ), intent(in) :: cmsg

  character(len=1, kind=c_char), dimension(:), pointer :: msg_fptr
  character(len=128) :: msg
  integer :: i, n

  call sofi_get_err_msg( int(cerr), msg )
  n = len_trim(msg)
  call c_f_pointer( cmsg, msg_fptr, [n+1] )
  do i = 1, n
     msg_fptr(i) = msg(i:i)
  end do
  msg_fptr(n+1) = c_null_char

end subroutine libira_get_err_msg
