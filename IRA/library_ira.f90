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
!! This file defines the wrappers to IRA routines, which are part of the
!! shared library shlib_ira.so.
!! The routines here are defined with "bind(C)" and "use iso_c_binding", and
!! they effectively behave as C-code. The arrays passed to these routines are
!! assumed to be already allocated by the caller, and are assumed to have
!! C-shape. Take care of transposed matrices, starting indices (1 or 0), etc.
!! Likewise, the output arrays are assumed to be already allocated by the caller,
!! and it's up to the caller to assure the correct shape and indices.
!! The routines here receive C-pointers to the memory where output should be
!! stored. Therefore, the output data appears as "intent(in)".
!!==============================================================================

subroutine lib_cshda( nat1, typ1, coords1, nat2, typ2, coords2, thr, found, dists )bind(C)
  !! wrapper to the cshda routine from cshda.f90
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat1
  type( c_ptr ), value, intent(in) :: typ1
  type( c_ptr ), value, intent(in) :: coords1
  integer(c_int), value, intent(in) :: nat2
  type( c_ptr ), value, intent(in) :: typ2
  type( c_ptr ), value, intent(in) :: coords2
  real( c_double ), value, intent(in) :: thr
  !!
  type( c_ptr ), intent(in) :: found
  type( c_ptr ), intent(in) :: dists

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_found
  type( c_ptr ), dimension(:), pointer :: p_atm1, p_atm2
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:), pointer :: p_dists

  !! local f memory
  real(c_double) :: some_thr


  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_atm1, [nat1] )
  call c_f_pointer( coords2, p_atm2, [nat2] )
  call c_f_pointer( p_atm1(1), p_coords1, [3,nat1] )
  call c_f_pointer( p_atm2(1), p_coords2, [3,nat2] )

  call c_f_pointer( found, p_found, [nat2] )
  call c_f_pointer( dists, p_dists, [nat2] )

  some_thr=thr

  call cshda( nat1, p_typ1, p_coords1, nat2, p_typ2, p_coords2, some_thr, &
       p_found, p_dists )

  !! output C-style indices, starting at 0
  p_found = p_found - 1

end subroutine lib_cshda



subroutine lib_cshda_pbc( nat1, typ1, coords1, nat2, typ2, coords2, lat, thr, found, dists )bind(C)
  !! wrapper to the cshda_pbc routine from cshda.f90
  use iso_c_binding
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
  type( c_ptr ), dimension(:), pointer :: p_atm1, p_atm2
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:), pointer :: p_dists, p_lat

  !! local f memory
  real(c_double), dimension(3,3) :: f_lat
  real(c_double) :: some_thr


  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_atm1, [nat1] )
  call c_f_pointer( coords2, p_atm2, [nat2] )
  call c_f_pointer( p_atm1(1), p_coords1, [3,nat1] )
  call c_f_pointer( p_atm2(1), p_coords2, [3,nat2] )
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

end subroutine lib_cshda_pbc


subroutine lib_ira_unify( nat1, typ1, coords1, candidate1, &
     nat2, typ2, coords2, candidate2, &
     kmax_factor, rotation, translation, permutation, hd )bind(C)
  !! wrapper to ira_unify routine from ira_routines.f90
  !!
  !! Warning:
  !! the indices in candidate1, and candidate2 need to be F-style (start by 1) on input!
  !!
  !! The result can be applied to struc2, as:
  !!
  !!   idx = permutation(i)
  !!   coords2(:,i) = matmul( rotation, coords2(:,idx) ) + translation
  !!
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat1
  type( c_ptr ), value, intent(in) :: typ1
  type( c_ptr ), value, intent(in) :: coords1
  type( c_ptr ), value, intent(in) :: candidate1
  integer(c_int), value, intent(in) :: nat2
  type( c_ptr ), value, intent(in) :: typ2
  type( c_ptr ), value, intent(in) :: coords2
  type( c_ptr ), value, intent(in) :: candidate2
  real( c_double ), value, intent(in) :: kmax_factor
  !!
  type( c_ptr ), intent(in) :: rotation
  type( c_ptr ), intent(in) :: translation
  type( c_ptr ), intent(in) :: permutation
  real(c_double), intent(out) :: hd

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_c1, p_c2, p_perm
  type( c_ptr ), dimension(:), pointer :: p_atm1, p_atm2
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:,:), pointer :: p_matrix
  real(c_double), dimension(:), pointer :: p_tr

  integer :: i, j

  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_atm1, [nat1] )
  call c_f_pointer( coords2, p_atm2, [nat2] )
  call c_f_pointer( p_atm1(1), p_coords1, [3,nat1] )
  call c_f_pointer( p_atm2(1), p_coords2, [3,nat2] )
  call c_f_pointer( candidate1, p_c1, [nat1] )
  call c_f_pointer( candidate2, p_c2, [nat2] )

  call c_f_pointer( rotation, p_matrix, [3,3] )
  call c_f_pointer( translation, p_tr, [3] )
  call c_f_pointer( permutation, p_perm, [nat2] )

  call ira_unify( nat1, p_typ1, p_coords1, p_c1, &
       nat2, p_typ2, p_coords2, p_c2, &
       kmax_factor, p_matrix, p_tr, p_perm, hd )

  ! write(*,*) 'Fortran matrix:'
  ! write(*,*) p_matrix(1,:)
  ! write(*,*) p_matrix(2,:)
  ! write(*,*) p_matrix(3,:)

  ! write(*,*) 'tr'
  ! write(*,*) p_tr

  ! write(*,*) nat1
  ! write(*,*)
  ! do i = 1, nat1
  !    write(*,*) p_typ1(i), p_coords1(:,i)
  ! end do
  ! write(*,*) nat2
  ! write(*,*)
  ! do i = 1, nat2
  !    j = p_perm(i)
  !    write(*,*) p_typ2(j), matmul( p_matrix, p_coords2(:,j) ) + p_tr
  ! end do

  !! return C-style data
  p_matrix = transpose( p_matrix )
  p_perm(:) = p_perm(:) - 1

end subroutine lib_ira_unify
