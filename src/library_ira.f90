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

!> @defgroup info_libira
!> @{

!!===============================================================================
!! @note
!! This file defines the wrappers to IRA routines, which are part of the
!! shared library libira.so.
!! The routines here are defined with "bind(C)" and "use iso_c_binding", and
!! they effectively behave as C-code. The arrays passed to these routines are
!! assumed to be already allocated by the caller, and are assumed to have
!! C-shape. Take care of transposed matrices, starting indices (1 or 0), etc.
!! Likewise, the output arrays are assumed to be already allocated by the caller,
!! and it's up to the caller to assure the correct shape and indices.
!! The routines here receive C-pointers to the memory where output should be
!! stored. Therefore, the output data appears as "intent(in)".
!!===============================================================================
!> @}

!> @brief wrapper to the cshda routine from cshda.f90
!!
!! @details
!! All parameters are as intent(in). The output is written into arrays
!! ``found`` and ``dists``, which are assumed to be allocated by the caller
!! to their correct size. ``size(found) = size(dists) = [nat2]``
!!
!! @param[in]  nat1    :: number of atoms in structure 1
!! @param[in]  typ1(nat1)    :: atomic types of structure 1
!! @param[in]  coords1(3,nat1) :: atomic positions of structure 1
!! @param[in]  nat2    :: number of atoms in structure 2
!! @param[in]  typ2(nat2)    :: atomic types of structure 2
!! @param[in]  coords2(3,nat2) :: atomic positions of structure 2
!! @param[in]  thr     :: threshold for the Hausdorff distance, used for early exit;
!! @param[in]  found(nat2)   :: list of assigned atoms of conf 2 to conf 1:
!!                              e.g. found(3) = 9 means atom 3 from conf 1 is assigned
!!                              to atom 9 in conf 2;
!! @param[in]  dists(nat2)   :: distances from atom i in conf 1 to atom found(i) in conf 2;
!! @return found, dists
!!
!!
!! C-header:
!! ~~~~~~~~~~~~~~~{.c}
!! void libira_cshda( int nat1, int *typ1, double *coords1, \
!!                 int nat2, int *typ2, double *coords2, \
!!                 double thr, int **found, double **dists);
!! ~~~~~~~~~~~~~~~
!!
subroutine libira_cshda( nat1, typ1, coords1, nat2, typ2, coords2, thr, found, dists )&
     bind(C, name="libira_cshda")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int),   value, intent(in) :: nat1
  type( c_ptr ),    value, intent(in) :: typ1
  type( c_ptr ),    value, intent(in) :: coords1
  integer(c_int),   value, intent(in) :: nat2
  type( c_ptr ),    value, intent(in) :: typ2
  type( c_ptr ),    value, intent(in) :: coords2
  real( c_double ), value, intent(in) :: thr
  !!
  type( c_ptr ), intent(in) :: found
  type( c_ptr ), intent(in) :: dists

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_found
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:), pointer :: p_dists

  !! local f memory
  real(c_double) :: some_thr

  integer :: i

  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )

  call c_f_pointer( coords1, p_coords1, [3,nat1] )
  call c_f_pointer( coords2, p_coords2, [3,nat2] )

  call c_f_pointer( found, p_found, [nat2] )
  call c_f_pointer( dists, p_dists, [nat2] )

  some_thr=thr

  call cshda( nat1, p_typ1, p_coords1, nat2, p_typ2, p_coords2, some_thr, &
       p_found, p_dists )

  !! output C-style indices, starting at 0
  p_found = p_found - 1

end subroutine libira_cshda


!> @brief wrapper to the cshda_pbc routine from cshda.f90
!!
!! @details
!! All parameters are as intent(in). The output is written into arrays
!! ``found`` and ``dists``, which are assumed to be allocated by the caller
!! to their correct size. ``size(found) = size(dists) = [nat2]``
!!
!! @param[in]  nat1    :: number of atoms in structure 1
!! @param[in]  typ1(nat1)    :: atomic types of structure 1
!! @param[in]  coords1(3,nat1) :: atomic positions of structure 1
!! @param[in]  nat2    :: number of atoms in structure 2
!! @param[in]  typ2(nat2)    :: atomic types of structure 2
!! @param[in]  coords2(3,nat2) :: atomic positions of structure 2
!! @param[in]  lat2(3,3)    :: lattice vectors of conf 2 in C order;
!! @param[in]  thr     :: threshold for the Hausdorff distance, used for early exit;
!! @param[in]  found(nat2)   :: list of assigned atoms of conf 2 to conf 1:
!!                        e.g. found(3) = 9 means atom 3 from conf 1 is assigned
!!                        to atom 9 in conf 2. Indices in C order (start at 0);
!! @param[in]  dists(nat2)   :: distances from atom i in conf 1 to atom found(i) in conf 2;
!! @return found, dists
!!
!!
!! C-header:
!! ~~~~~~~~~~~~~~~{.c}
!! void libira_cshda_pbc( int nat1, int *typ1, double *coords1, \
!!                     int nat2, int *typ2, double *coords2, double * lat2, \
!!                     double thr, int **found, double **dists);
!! ~~~~~~~~~~~~~~~
!!
subroutine libira_cshda_pbc( nat1, typ1, coords1, nat2, typ2, coords2, lat, thr, found, dists )&
     bind(C, name="libira_cshda_pbc")
  !! wrapper to the cshda_pbc routine from cshda.f90
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat1
  type( c_ptr ), value, intent(in) :: typ1
  type( c_ptr ), value, intent(in) :: coords1
  integer(c_int), value, intent(in) :: nat2
  type( c_ptr ), value, intent(in) :: typ2
  type( c_ptr ), value, intent(in) :: coords2
  type( c_ptr ), value, intent(in) :: lat
  real( c_double ), value, intent(in) :: thr
  !!
  type( c_ptr ), intent(in) :: found
  type( c_ptr ), intent(in) :: dists

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_found
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:), pointer :: p_dists, p_lat

  !! local f memory
  real(c_double), dimension(3,3) :: f_lat
  real(c_double) :: some_thr


  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_coords1, [3,nat1] )
  call c_f_pointer( coords2, p_coords2, [3,nat2] )
  call c_f_pointer( lat, p_lat, [9] )
  f_lat = reshape( p_lat, [3,3] )
  !! receive C-style array, do transpose
  f_lat = transpose( f_lat )

  ! write(*,*) "Received lattice vectors:"
  ! write(*,*) "vec 1:", f_lat(1,:)
  ! write(*,*) "vec 2:", f_lat(2,:)
  ! write(*,*) "vec 3:", f_lat(3,:)

  call c_f_pointer( found, p_found, [nat2] )
  call c_f_pointer( dists, p_dists, [nat2] )

  some_thr=thr

  call cshda_pbc( nat1, p_typ1, p_coords1, nat2, p_typ2, p_coords2, f_lat, &
       some_thr, p_found, p_dists )

  !! output C-style array, strating indices at 0
  p_found = p_found - 1

end subroutine libira_cshda_pbc

!> @brief the IRA matching procedure
!!
!! @details
!! wrapper call to match two structures. This includes call to ira_unify to get apx,
!! and then call to svdrot_m to obtain final match. Routines from ira_routines.f90
!!
!! @note
!!  Warning:
!!  the indices in candidate1, and candidate2 need to be F-style (start by 1) on input!
!!
!! The result can be applied to struc2, as:
!! ~~~~~~~~~~~~~~~{.f90}
!!    idx = permutation(i)
!!    coords2(:,i) = matmul( rotation, coords2(:,idx) ) + translation
!! ~~~~~~~~~~~~~~~
!!
!! @param[in]  nat1    :: number of atoms in structure 1
!! @param[in]  typ1(nat1)    :: atomic types of structure 1
!! @param[in]  coords1(3,nat1) :: atomic positions of structure 1
!! @param[in]  candidate1(nat1) :: list of candidate central atoms in structure 1
!! @param[in]  nat2    :: number of atoms in structure 2
!! @param[in]  typ2(nat2)    :: atomic types of structure 2
!! @param[in]  coords2(3,nat2) :: atomic positions of structure 2
!! @param[in]  candidate2(nat2) :: list of candidate central atoms in structure 2
!! @param[in]  kmax_factor :: the factor to multiply kmax (should be > 1.0)
!! @param[in]  rotation(3,3) :: the 3x3 rotation matrix in C order
!! @param[in]  translation(3) :: the 3D translation vector
!! @param[in]  permutation(nat2) :: the atomic permutations in C order (start at 0)
!! @param[out]  hd :: final value of the Hausdorff distance
!! @param[out]  cerr :: error value (negative on error, zero otherwise)
!! @returns rotation, translation, permutation, hd, cerr
!!
!! C-header:
!! ~~~~~~~~~~~~~~~{.c}
!! void libira_match(int nat1, int *typ1, double *coords1, int *cand1,\
!!                int nat2, int *typ2, double *coords2, int *cand2, \
!!                double km_factor, double **rmat, double **tr, \
!!                int **perm, double *hd, int *cerr);
!! ~~~~~~~~~~~~~~~
subroutine libira_match( nat1, typ1, coords1, candidate1, &
     nat2, typ2, coords2, candidate2, &
     kmax_factor, rotation, translation, permutation, hd, cerr ) bind(C, name="libira_match")
  use, intrinsic :: iso_c_binding
  use err_module, only: get_err_msg
  implicit none
  integer(c_int), value, intent(in) :: nat1
  type( c_ptr ), value,  intent(in) :: typ1
  type( c_ptr ), value,  intent(in) :: coords1
  type( c_ptr ), value,  intent(in) :: candidate1
  integer(c_int), value, intent(in) :: nat2
  type( c_ptr ), value,  intent(in) :: typ2
  type( c_ptr ), value,  intent(in) :: coords2
  type( c_ptr ), value,  intent(in) :: candidate2
  real( c_double ), value, intent(in) :: kmax_factor
  !!
  type( c_ptr ), intent(in) :: rotation
  type( c_ptr ), intent(in) :: translation
  type( c_ptr ), intent(in) :: permutation
  real(c_double), intent(out) :: hd
  integer( c_int ), intent(out) :: cerr

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_c1, p_c2, p_perm
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:,:), pointer :: p_matrix
  real(c_double), dimension(:), pointer :: p_tr

  integer :: i, ierr
  real( c_double ), dimension(3,nat2) :: fcoords2
  integer(c_int), dimension(nat2) :: ftyp2
  real( c_double ), dimension(3,3) :: srot
  real( c_double ), dimension(3) :: str, rdum
  real( c_double ) :: pp

  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_coords1, [3,nat1] )
  call c_f_pointer( coords2, p_coords2, [3,nat2] )
  call c_f_pointer( candidate1, p_c1, [nat1] )
  call c_f_pointer( candidate2, p_c2, [nat2] )

  call c_f_pointer( rotation, p_matrix, [3,3] )
  call c_f_pointer( translation, p_tr, [3] )
  call c_f_pointer( permutation, p_perm, [nat2] )

  !! get apx
  call ira_unify( nat1, p_typ1, p_coords1, p_c1, &
                  nat2, p_typ2, p_coords2, p_c2, &
                  kmax_factor, p_matrix, p_tr, p_perm, hd, ierr )
  cerr = int( ierr, c_int )
  ! write(*,*) "HD after unify",hd
  if( ierr /= 0 ) then
     write(*,*) "ERROR in libira_match"
     write(*,*) get_err_msg( ierr )
     return
  end if


  !! transform
  ftyp2(:) = p_typ2(p_perm(:))
  fcoords2(:,:) = p_coords2(:,p_perm(:))
  do i = 1, nat2
     fcoords2(:,i) = matmul( p_matrix, fcoords2(:,i)) + p_tr
  end do

  !! call svd
  call svdrot_m( nat1, p_typ1, p_coords1, &
       nat1, ftyp2(1:nat1), fcoords2(:,1:nat1), &
       srot, str, ierr )
  if( ierr /= 0 ) then
     return
  end if


  !! apply svd
  do i = 1, nat2
     fcoords2(:,i) = matmul( srot, fcoords2(:,i) ) + str
  end do

  !! measure dH
  pp=0.0_c_double
  do i = 1, nat1
     rdum = p_coords1(:,i) - fcoords2(:,i)
     pp = max( pp, norm2(rdum) )
  end do
  hd=pp

  !! put together
  p_matrix = matmul( srot, p_matrix )
  p_tr = matmul(srot, p_tr) + str

  !! return C-style data
  p_matrix = transpose( p_matrix )
  p_perm(:) = p_perm(:) - 1

end subroutine libira_match


!> @details
!! Get the IRA version string, and date.
!! This is a wrapper to get_version from version.f90
!!
!! C header:
!!~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
!! void libira_get_version( char *string, int *date );
!!~~~~~~~~~~~~~~~~~~~~~~~~~
!!
!! @param[out] cstring(10) :: version string
!! @param[out] cdate :: version date, format YYYYmm
!! @returns cstring, cdate
!!
subroutine libira_get_version( cstring, cdate )bind(C, name="libira_get_version")
  use, intrinsic :: iso_c_binding, only: c_int, c_null_char, c_char
  implicit none
  character(len=1, kind=c_char), dimension(6) :: cstring
  integer( c_int ) :: cdate

  character(5) :: fstring
  integer :: fdate, i, n

  call ira_get_version( fstring, fdate )
  cdate = int( fdate, c_int )
  n = len_trim(fstring)
  if( n .gt. 5 ) write(*,*) "WARNING: IRA version string seems long, check!"
  do i = 1, n
     cstring(i) = fstring(i:i)
  end do
  cstring(n+1) = c_null_char

end subroutine libira_get_version
