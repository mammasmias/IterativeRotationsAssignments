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

subroutine lib_compute_all( nat, typ, coords, sym_thr, &
                            nmat, mat_list, perm_list, &
                            op_list, n_list, p_list, &
                            ax_list, angle_list, dmax_list, pg ) bind(C)
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  integer( c_int ), value, intent(in) :: nat
  type( c_ptr ), value, intent(in) :: typ
  type( c_ptr ), value, intent(in) :: coords
  real( c_double ), value, intent(in) :: sym_thr

  integer( c_int ), intent(out) :: nmat
  type( c_ptr ), intent(in) :: mat_list
  type( c_ptr ), intent(in) :: perm_list
  type( c_ptr ), intent(in) :: op_list
  type( c_ptr ), intent(in) :: n_list
  type( c_ptr ), intent(in) :: p_list
  type( c_ptr ), intent(in) :: ax_list
  type( c_ptr ), intent(in) :: angle_list
  type( c_ptr ), intent(in) :: dmax_list
  type( c_ptr ), intent(in) :: pg

  !! pointers for c input arrays
  integer(c_int), dimension(:), pointer :: ptyp
  real( c_double), dimension(:,:), pointer :: pcoords
  type( c_ptr ), dimension(:), pointer :: patom
  real( c_double ), dimension(:,:,:), pointer :: pmat_list
  integer( c_int ), dimension(:,:), pointer :: pperm_list
  integer( c_int ), dimension(:), pointer :: pn_list
  integer( c_int ), dimension(:), pointer :: pp_list
  real( c_double ), dimension(:,:), pointer :: pax_list
  real( c_double ), dimension(:), pointer :: pangle_list
  real( c_double ), dimension(:), pointer :: pdmax_list
  character(len=1, kind=c_char), dimension(:), pointer :: pg_char
  character(len=1, kind=c_char), dimension(:), pointer :: op_char

  !! some f-defined memory
  integer( c_int ) :: nb
  character(len=10) :: f_pg
  character(len=2), dimension(nmax) :: fop_list
  character(len=2) :: this_op

  integer :: i, m, n, lenc

  ! write(*,*) "entering lib"
  !! local arrays for computation

  !!
  !! receive input
  !!
  call c_f_pointer( typ, ptyp, [nat] )
  call c_f_pointer( coords, patom, [nat] )
  call c_f_pointer( patom(1), pcoords, [3,nat] )
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

  !!
  !! compute SOFI
  !!
  call sofi_compute_all( nat, ptyp, pcoords, sym_thr, &
       nb, pmat_list, pperm_list, &
       fop_list, pn_list, pp_list, &
       pax_list, pangle_list, pdmax_list, f_pg )

  !! set pg string
  lenc=len_trim(f_pg)
  do i=1, lenc
     pg_char(i) = f_pg(i:i)
  end do
  pg_char(lenc+1)=c_null_char

  !! compact the op_list into single string to pass to C
  lenc=2*nmax
  call c_f_pointer( op_list, op_char, [lenc+1] )
  n=1
  do i = 1, nb
     this_op=fop_list(i)
     do m = 1, 2
        op_char(n) = this_op(m:m)
        n = n + 1
     end do
  end do
  op_char(n)=c_null_char

  !! set output value for nmat
  nmat = nb

  !! transpose all output matrices to C-order
  do i = 1, nb
     pmat_list(:,:,i) = transpose( pmat_list(:,:,i))
  end do

  !! return C-style indices
  pperm_list = pperm_list-1

  ! write(*,*) "exiting lib"
end subroutine lib_compute_all


subroutine lib_get_symm_ops(nat, typ, coords, symm_thr, n_sym, sym_list )bind(C)
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  !! "input" structure
  integer( c_int ), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  real( c_double ), value, intent(in) :: symm_thr
  !! "output"
  integer( c_int ), intent(out) :: n_sym
  type( c_ptr ),    intent(in) :: sym_list

  !! F pointers
  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  type( c_ptr ), dimension(:), pointer :: ptr_atom
  real( c_double), dimension(:,:,:), pointer :: ptr_op
  !! memory in f
  integer(c_int) :: n, i

  ! write(*,*) 'enter symmop'
  !!
  !! receive input
  !!
  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_atom, [nat] )
  call c_f_pointer( ptr_atom(1), ptr_coords, [3,nat] )

  call c_f_pointer( sym_list, ptr_op, [3,3,nmax] )

  call sofi_get_symmops( nat, ptr_typ, ptr_coords, symm_thr, n, ptr_op )

  !! set output data
  n_sym=n

  !! transpose output matrices to C-order
  do i = 1, n
     ptr_op(:,:,i) = transpose( ptr_op(:,:,i) )
  end do


end subroutine lib_get_symm_ops


subroutine lib_get_pg( nbas, cptr_op_list, ppg )bind(C)
  use iso_c_binding
  implicit none
  integer( c_int ), value, intent(in) :: nbas
  type( c_ptr ), value, intent(in) :: cptr_op_list
  !! c ptr to write output
  type( c_ptr ), intent(in) :: ppg

  !! f pointers
  type( c_ptr ), dimension(:), pointer :: p_lvl1
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real(c_double), dimension(:,:,:), pointer :: op_list
  ! character(len=10, kind=c_char), pointer :: str
  character(len=1, kind=c_char), dimension(:), pointer :: pg_char
  !! local
  integer :: i, n
  character(len=10) :: pg
  logical :: verb

  !! receive input
  call c_f_pointer( cptr_op_list, p_lvl1, [nbas] )
  call c_f_pointer( p_lvl1(1), p_lvl2, [9,nbas] )


  !! put the array into proper fortran shape
  allocate( op_list(1:3,1:3,1:nbas))
  do n = 1, nbas
     op_list(:,:,n)=reshape( p_lvl2(:,n),[3,3])
     !! transpose because it is C-style on input
     op_list(:,:,n) = transpose( op_list(:,:,n) )
  end do

  call c_f_pointer( ppg, pg_char, [11] )

  verb = .false.
  call sofi_get_pg( nbas, op_list, pg, verb )

  n = len_trim(pg)
  do i = 1, n
     pg_char(i)=pg(i:i)
  end do
  pg_char(n+1)=c_null_char

end subroutine lib_get_pg


subroutine lib_unique_ax_angle( n_mat, cptr_mat_list, op_out, ax_out, angle_out )bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: n_mat
  type( c_ptr ), value, intent(in) :: cptr_mat_list

  type( c_ptr ), intent(in) :: op_out
  type( c_ptr ), intent(in) :: ax_out
  type( c_ptr ), intent(in) :: angle_out

  type( c_ptr ), dimension(:), pointer :: p_lvl1
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real(c_double), dimension(:,:,:), pointer :: mat_list
  real( c_double ), dimension(:,:), pointer :: ptr_ax
  real( c_double ), dimension(:), pointer :: ptr_angle
  character(len=1, kind=c_char), dimension(:), pointer :: op_char

  integer :: i, n, m
  character(len=2), dimension(n_mat) :: fop_list
  character(len=2) :: this_op

  !! receive input
  call c_f_pointer( cptr_mat_list, p_lvl1, [n_mat] )
  call c_f_pointer( p_lvl1(1), p_lvl2, [9,n_mat] )
  call c_f_pointer( ax_out, ptr_ax, [3,n_mat] )
  call c_f_pointer( angle_out, ptr_angle, [n_mat] )

  !! put arrays into proder fortran order
  allocate( mat_list(1:3,1:3,1:n_mat))
  do n = 1, n_mat
     mat_list(:,:,n)=reshape( p_lvl2(:,n),[3,3])
     !! transpose becuse it's C-order on input
     mat_list(:,:,n) = transpose(mat_list(:,:,n) )
  end do

  call sofi_unique_ax_angle( n_mat, mat_list, fop_list, ptr_ax, ptr_angle )

  !! concatenate op strings into single string len(2*nmat+1)
  n=2*n_mat
  call c_f_pointer( op_out, op_char, [n+1] )
  n=1
  do i = 1, n_mat
     this_op = fop_list(i)
     do m = 1, 2
        op_char(n) = this_op(m:m)
        n = n + 1
     end do
  end do
  op_char(n) = c_null_char

end subroutine lib_unique_ax_angle


subroutine lib_analmat( c_rmat, c_op, n, p, c_ax, angle )bind(C)
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

  !! f pointers
  type( c_ptr ), dimension(:), pointer :: ax1
  real(c_double), dimension(:,:), pointer :: rmat
  real( c_double), dimension(:), pointer :: ax
  character(len=1, kind=c_char), dimension(:), pointer :: op_fptr

  !! local mem
  integer :: i, m
  character(len=2) :: op

  !! connect c to f
  call c_f_pointer( c_rmat, ax1, [3] )
  call c_f_pointer( ax1(1), rmat, [3,3] )
  rmat = transpose(rmat)
  call c_f_pointer( c_op, op_fptr, [3] )
  call c_f_pointer( c_ax, ax, [3] )

  call sofi_analmat( rmat, op, n, p, ax, angle )
  m=len_trim(op)
  do i = 1, m
     op_fptr(i) = op(i:i)
  end do
  op_fptr(m+1)=c_null_char

end subroutine lib_analmat


subroutine lib_ext_bfield( n_mat, cop_list, cb_field, n_out, cop_out )bind(C)
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
  type( c_ptr ), dimension(:), pointer :: p_lvl1
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double ), dimension(:), pointer :: b_field

  !! f memory
  real(c_double), dimension(:,:,:), pointer :: pmat_list
  real(c_double), dimension(3,3,n_mat) :: mat_list
  integer :: n
  integer :: m, i

  m=n_mat

  !! connect c to f
  call c_f_pointer( cop_list, p_lvl1, [n_mat] )
  call c_f_pointer( p_lvl1(1), p_lvl2, [9,n_mat] )

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



end subroutine lib_ext_bfield


subroutine lib_get_perm( nat, typ, coords, nbas, bas_list, perm_list, dmax_list)bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  integer( c_int ), value, intent(in) :: nbas
  type( c_ptr ), value,  intent(in) :: bas_list
  !! "output"
  type( c_ptr ), intent(in) :: perm_list
  type( c_ptr ), intent(in) :: dmax_list

  !! F pointers
  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  type( c_ptr ), dimension(:), pointer :: ptr_atom
  type( c_ptr ), dimension(:), pointer :: p_lvl1
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double), dimension(:,:,:), pointer :: ptr_op
  integer( c_int ), dimension(:,:), pointer :: pperm_list
  real( c_double ), dimension(:), pointer :: pdmax_list

  integer :: i

  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_atom, [nat] )
  call c_f_pointer( ptr_atom(1), ptr_coords, [3,nat] )

  call c_f_pointer( perm_list, pperm_list, [nat,nbas] )
  call c_f_pointer( dmax_list, pdmax_list, [nbas] )

  call c_f_pointer( bas_list, p_lvl1, [nbas])
  call c_f_pointer( p_lvl1(1), p_lvl2, [9,nbas] )

  !! receive input array and reshape into fortran order
  allocate( ptr_op(1:3,1:3,1:nbas))
  do i = 1, nbas
     ptr_op(:,:,i) = reshape( p_lvl2(:,i),[3,3])
     !! transpose because it is C-order on input
     ptr_op(:,:,i) = transpose( ptr_op(:,:,i) )
  end do

  call sofi_get_perm( nat, ptr_typ, ptr_coords, nbas, ptr_op, pperm_list, pdmax_list )

  !! output C-style: start at 0 indices
  pperm_list = pperm_list - 1

end subroutine lib_get_perm


subroutine lib_get_combos( nat, typ, coords, nbas_in, bas_in, nbas_out, bas_out )bind(C)
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  integer(c_int), value, intent(in) :: nat
  type( c_ptr ),    value, intent(in) :: typ
  type( c_ptr ),    value, intent(in) :: coords
  integer( c_int ), value, intent(in) :: nbas_in
  type( c_ptr ), value,  intent(in) :: bas_in
  integer( c_int ), intent(out) :: nbas_out
  type( c_ptr ), intent(in) :: bas_out


  integer(c_int), dimension(:), pointer :: ptr_typ
  real(c_double), dimension(:,:), pointer :: ptr_coords
  type( c_ptr ), dimension(:), pointer :: ptr_atom
  type( c_ptr ), dimension(:), pointer :: p_lvl1
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double), dimension(:,:,:), pointer :: pcombo_out

  real(c_double), dimension(3,3,nmax) :: bas_list
  integer :: i, m

  call c_f_pointer( typ, ptr_typ, [nat] )
  call c_f_pointer( coords, ptr_atom, [nat] )
  call c_f_pointer( ptr_atom(1), ptr_coords, [3,nat] )

  call c_f_pointer( bas_in, p_lvl1, [nbas_in])
  call c_f_pointer( p_lvl1(1), p_lvl2, [9,nbas_in] )

  do i = 1, nbas_in
     bas_list(:,:,i) = reshape( p_lvl2(:,i),[3,3])
     bas_list(:,:,i) = transpose(bas_list(:,:,i))
  end do

  m = nbas_in

  call sofi_get_combos( nat, ptr_typ, ptr_coords, m, bas_list )

  nbas_out = m
  call c_f_pointer( bas_out, pcombo_out, [3,3,m] )
  pcombo_out = bas_list(1:3,1:3,1:m)


end subroutine lib_get_combos


subroutine lib_try_mat( nat, typ, coords, rmat, dh, perm )bind(C)
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
  type( c_ptr ), dimension(:), pointer :: ptr_atom
  type( c_ptr ), dimension(:), pointer :: p_lvl1
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
  call c_f_pointer( coords, ptr_atom, [nat] )
  call c_f_pointer( ptr_atom(1), ptr_coords, [3,nat] )

  call c_f_pointer( rmat, p_lvl1, [3] )
  call c_f_pointer( p_lvl1(1), p_lvl2, [3,3] )
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


end subroutine lib_try_mat


subroutine lib_construct_operation( op, axis, angle, matrix )bind(C)
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

  real( c_double ), dimension(:), pointer :: ptr_ax
  real(c_double), dimension(:,:), pointer :: ptr_matrix

  real( c_double ), dimension(3,3) :: f_matrix
  character(len=1,kind=c_char), dimension(:), pointer :: ptr_op
  character(len=2) :: f_op
  integer( c_int ) :: i, n

  n=c_strlen(op)
  call c_f_pointer( op, ptr_op, [n] )
  f_op=""
  do i = 1, n
     f_op(i:i)=ptr_op(i)
  end do

  f_op=trim(adjustl(f_op))

  call c_f_pointer( axis, ptr_ax, [3] )
  call c_f_pointer( matrix, ptr_matrix, [3,3] )

  call sofi_construct_operation( f_op, ptr_ax, angle, f_matrix )

  !! return C-order matrix
  ptr_matrix = transpose(f_matrix)

end subroutine lib_construct_operation


subroutine lib_mat_combos( nbas_in, bas_in, nbas_out, bas_out )bind(C)
  use iso_c_binding
  use sofi_tools, only: nmax
  implicit none
  integer( c_int ), value, intent(in) :: nbas_in
  type( c_ptr ), value,  intent(in) :: bas_in
  integer( c_int ), intent(out) :: nbas_out
  type( c_ptr ), intent(in) :: bas_out


  type( c_ptr ), dimension(:), pointer :: p_lvl1
  real( c_double ), dimension(:,:), pointer :: p_lvl2
  real( c_double), dimension(:,:,:), pointer :: pcombo_out

  real(c_double), dimension(3,3,nmax) :: bas_list
  real(c_double), dimension(3,3,nbas_in) :: bas_inp
  integer :: i, m

  call c_f_pointer( bas_in, p_lvl1, [nbas_in])
  call c_f_pointer( p_lvl1(1), p_lvl2, [9,nbas_in] )

  do i = 1, nbas_in
     bas_list(:,:,i) = reshape( p_lvl2(:,i),[3,3])
     bas_list(:,:,i) = transpose(bas_list(:,:,i))
  end do

  m = nbas_in

  bas_inp(:,:,:) = bas_list(:,:,1:nbas_in)

  call sofi_mat_combos( nbas_in, bas_inp, m, bas_list )

  nbas_out = m
  call c_f_pointer( bas_out, pcombo_out, [3,3,m] )
  pcombo_out = bas_list(1:3,1:3,1:m)


end subroutine lib_mat_combos
