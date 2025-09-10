!!
!! Copyright (C) 2021, MAMMASMIAS Consortium
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

!> @brief determinant of a 3x3 matrix a
subroutine determinant3x3(a, d)
  use ira_precision
  implicit none
  real(rp), dimension(3, 3), intent(in) :: a
  real(rp), intent(out) :: d

  d = a(1, 1)*a(2, 2)*a(3, 3) - a(1, 1)*a(2, 3)*a(3, 2) + &
       a(1, 2)*a(2, 3)*a(3, 1) - a(2, 2)*a(3, 1)*a(1, 3) + &
       a(2, 1)*a(3, 2)*a(1, 3) - a(3, 3)*a(1, 2)*a(2, 1)

end subroutine determinant3x3

!> @brief Use mergesort to sort the input array by the chosen axis.
!! @details The routine mergesort is located in sorting_module.
subroutine sort(n, ndim, array, axis)
  use sorting_module, only: mergesort
  use ira_precision
  implicit none
  integer(ip), intent(in) :: n
  integer(ip), intent(in) :: ndim
  real(rp), dimension(ndim, n), intent(inout) :: array
  integer(ip), intent(in) :: axis

  real(rp), allocatable :: work(:, :)

  !! allocate working space for the sorting routine
  allocate (work(1:ndim, 1:(n + 1)/2), source=0.0_rp)

  call mergesort(array, work, axis)

  deallocate (work)

end subroutine sort


!> @brief routine that calls LAPACK svd routine
!!
!! @param[in] m     :: leading dimension of matrix a
!! @param[in] n     :: second dimension of matrix a
!! @param[in] a     :: matrix a
!! @param[out] s     :: diagonal matrix of singular values
!! @param[out] u     :: orthonormal matrix u
!! @param[out] v     :: orthonormal matrix v
!! @param[out] ierr  :: error code
subroutine svd(m, n, a, u, s, v, ierr)
  use ira_precision
  use err_module
  implicit none
  interface
     ! lapack
     subroutine dgesvd( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
       use, intrinsic :: iso_fortran_env, only: ddp => real64
       character(len=1), intent(in) :: JOBU, JOBVT
       integer, intent(in)   :: LDA, LDU, LDVT, LWORK, M, N
       integer, intent(out)   :: INFO
       real(ddp), intent(inout) :: A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * )
       real(ddp), intent(out) :: WORK( * )
     end subroutine dgesvd
  end interface
  integer(ip), intent(in) :: m
  integer(ip), intent(in) :: n
  real(rp), dimension(m, n), intent(in) :: a
  real(rp), dimension(m, n), intent(out) :: s
  real(rp), dimension(m, n), intent(out) :: u
  real(rp), dimension(m, n), intent(out) :: v
  integer(ip), intent(out) :: ierr

  real(rp), dimension(m, n) :: a_copy
  real(rp), dimension(min(m, n)) :: sdiag
  real(rp), allocatable :: work(:)
  integer(ip) :: lwork, i, info, lda, ldu, ldv
  character(len=3) :: jobu, jobv

  ierr = 0

  lwork = max(3*min(m, n) + max(m, n), 5*min(m, n))

  allocate (work(1:lwork))
  !
  !  Compute the eigenvalues and eigenvectors.
  !
  jobu = 'A'
  jobv = 'A'
  lda = m
  ldu = m
  ldv = n
  !
  !  The input matrix is destroyed by the routine.  Since we need to keep
  !  it around, we only pass a copy to the routine.
  !
  a_copy(1:m, 1:n) = a(1:m, 1:n)

  !!
  !! for the single precision real use this:
  ! call sgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
  !     lwork, info )
  !!
  !!
  !! for the double precision real use this:
  call dgesvd(jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
       lwork, info)

  if (info /= 0) then
     write (*, '(a)') ' '
     write (*, '(a)') '  The SVD could not be calculated.'
     write (*, '(a)') '  LAPACK routine DGESVD returned a nonzero'
     write (*, '(a,i8)') '  value of the error flag, INFO = ', info
     ierr = ERR_SVD
     return
  end if
  !
  !  Make the MxN matrix S from the diagonal values in SDIAG.
  !
  s(1:m, 1:n) = 0.0_rp
  do i = 1, min(m, n)
     s(i, i) = sdiag(i)
  end do
  !
  !  Transpose V.
  !
  !    v = transpose ( v )

  deallocate (work)

  return

end subroutine svd

! subroutine cross_prod( a, b, c )
!   !> @brief Cross product of two vectors
!   implicit none
!   real(rp), dimension(3), intent(in) :: a
!   real(rp), dimension(3), intent(in) :: b
!   real(rp), dimension(3), intent(out) :: c

!   c(1) = a(2)*b(3) - a(3)*b(2)
!   c(2) = a(3)*b(1) - a(1)*b(3)
!   c(3) = a(1)*b(2) - a(2)*b(1)

! end subroutine cross_prod

!> @details
!! permute a real 2D array into order, equivalent to:
!!
!!   array(:,:) = array(:, order(:) )
!!
subroutine permute_real_2d(n, m, array, order)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: n
  integer(ip), intent(in) :: m
  real(rp), dimension(m, n), intent(inout) :: array
  integer(ip), dimension(n), intent(in) :: order

  integer(ip) :: i
  real(rp), dimension(m, n) :: tmp

  !! tmp copy
  tmp(:, :) = array(:, :)

  !! permute
  do i = 1, n
     array(:, i) = tmp(:, order(i))
  end do

end subroutine permute_real_2d

!> @details
!! permute a real 2D array into inverse order,
!! equivalent to:
!!
!!   array(:, order(:) ) = array(:,:)
!!
subroutine permute_real_2d_back(n, m, array, order)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: n
  integer(ip), intent(in) :: m
  real(rp), dimension(m, n), intent(inout) :: array
  integer(ip), dimension(n), intent(in) :: order

  integer(ip) :: i
  real(rp), dimension(m, n) :: tmp

  !! tmp copy
  tmp(:, :) = array(:, :)

  !! permute
  do i = 1, n
     array(:, order(i)) = tmp(:, i)
  end do

end subroutine permute_real_2d_back

!> @details
!! permute an integer 1D array into order, equivalent to:
!!
!!    array(:) = array( order(:) )
!!
subroutine permute_int_1d(n, array, order)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: n
  integer(ip), dimension(n), intent(inout) :: array
  integer(ip), dimension(n), intent(in) :: order

  integer(ip) :: i
  integer(ip), dimension(n) :: tmp

  !! tmp copy
  tmp(:) = array(:)

  !! permute
  do i = 1, n
     array(i) = tmp(order(i))
  end do

end subroutine permute_int_1d

!> @details
!! permute an integer 1D array into inverse order,
!! equivalent to:
!!
!!    array( order(:) ) = array(:)
!!
subroutine permute_int_1d_back(n, array, order)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: n
  integer(ip), dimension(n), intent(inout) :: array
  integer(ip), dimension(n), intent(in) :: order

  integer(ip) :: i
  integer(ip), dimension(n) :: tmp

  !! tmp copy
  tmp(:) = array(:)

  !! permute
  do i = 1, n
     array(order(i)) = tmp(i)
  end do

end subroutine permute_int_1d_back

!> @brief Set orthonormal basis from input vectors, third is cross of
!! the first two, fail happens when input is collinear.
!!
!! basis on output contains basis vectors in rows.
subroutine set_orthonorm_bas(vec1, vec2, basis, fail)
  use ira_precision
  use sofi_tools, only: collinearity_thr
  implicit none

  real(rp), dimension(3), intent(in) :: vec1, vec2
  real(rp), dimension(3, 3), intent(out) :: basis
  logical, intent(out) :: fail

  real(rp) :: prod
  real(rp) :: small_size_thr
  real(rp) :: norm_v

  !! threshold for too small vectors
  !! (for numerical reasons, don't want to normalize a too small vector)
  small_size_thr = 1e-1_rp

  fail = .true.
  basis(:, :) = 0.0_rp
  ! write(*,*) vec1
  ! write(*,*) vec2

  !! first vector, normalize
  norm_v = sqrt(dot_product(vec1, vec1))
  !! if too small, return with fail = .true.
  if (norm_v .lt. small_size_thr) return
  basis(1, :) = vec1(:)/norm_v

  !! second vector, normalize
  norm_v = sqrt(dot_product(vec2, vec2))
  !! if too small, return with fail = .true.
  if (norm_v .lt. small_size_thr) return
  basis(2, :) = vec2(:)/norm_v

  !! check projection
  prod = dot_product(basis(1, :), basis(2, :))

  !! if vectors are collinear, return with fail = .true.
  if (abs(prod) .gt. collinearity_thr) then
     fail = .true.
     ! write(*,*) 'failed collinearity',abs(prod)
     return
  end if

  !! take orthogonal component of second vector, normalize
  basis(2, :) = basis(2, :) - prod*basis(1, :)
  norm_v = sqrt(dot_product(basis(2, :), basis(2, :)))
  basis(2, :) = basis(2, :)/norm_v

  !! third vector is cross product of first two, normalize
  basis(3, 1) = basis(1, 2)*basis(2, 3) - basis(1, 3)*basis(2, 2)
  basis(3, 2) = basis(1, 3)*basis(2, 1) - basis(1, 1)*basis(2, 3)
  basis(3, 3) = basis(1, 1)*basis(2, 2) - basis(1, 2)*basis(2, 1)

  norm_v = sqrt(dot_product(basis(3, :), basis(3, :)))
  basis(3, :) = basis(3, :)/norm_v

  !!
  !! check for NaN. Can it happen?
  !!
  if (any(basis .ne. basis)) then
     fail = .true.
     ! write(*,*) 'fail NaN'
     ! write(*,*) basis(1,:)
     ! write(*,*) basis(2,:)
     ! write(*,*) basis(3,:)
     return
  end if

  !! if we come to this point, the basis is ok
  fail = .false.
  return
end subroutine set_orthonorm_bas

!> @details
!! Optimal rotation by SVD.
!!
!! This rouine allows rotation matrix to be a reflection.
!! no 'correction' when determinant is negative
!!
!! The output can be applied as:
!!
!!~~~~~~~~~~~~~{.f90}
!!     coords2(:,i) = matmul( rmat, coords2(:,i) ) + translate
!!~~~~~~~~~~~~~
!!
!! OR
!!
!!~~~~~~~~~~~~~~~{.f90}
!!     coords1(:,i) = matmul( transpose(rmat), coords1(:,i) ) - &
!!                                  matmul( transpose(rmat), translate )
!!~~~~~~~~~~~~~~~
!!
!!
!! @param[in] nat1    -> number of atoms in conf 1;
!! @param[in] typ1(nat1)    -> atomic types in conf 1;
!! @param[in] coords1_in(3,nat1) -> coordinates of conf 1;
!! @param[in] nat2    -> number of atoms in conf 2;
!! @param[in] typ2(nat2)    -> atomic types in conf 2;
!! @param[in] coords2_in(3,nat2) -> coordinates of conf 2;
!!
!! @param[out] rmat(3,3)    -> 3x3 optimal rotation matrix
!! @param[out] translate(3) -> optimal translation vector
!! @returns rmat, translate
!!
subroutine svdrot_m(nat1, typ1, coords1_in, &
     nat2, typ2, coords2_in, &
     rmat, translate, ierr)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  real(rp), dimension(3, 3), intent(out) :: rmat
  real(rp), dimension(3), intent(out) :: translate
  integer(ip), intent(out) :: ierr

  real(rp), dimension(3, nat1) :: coords1
  real(rp), dimension(3, nat2) :: coords2
  real(rp), dimension(3) :: gc1, gc2
  real(rp), dimension(3, 3) :: matrix, u, smat, vt
  integer(ip) :: i

  !! set initial copies
  coords1(:, :) = coords1_in(:, :)
  coords2(:, :) = coords2_in(:, :)

  !! set geo centers
  gc1(:) = 0.0_rp
  gc2(:) = 0.0_rp
  gc1(:) = sum(coords1(:, :), 2)/nat1
  gc2(:) = sum(coords2(:, :), 2)/nat2

  !! recenter
  do i = 1, nat1
     coords1(:, i) = coords1(:, i) - gc1
  end do
  do i = 1, nat2
     coords2(:, i) = coords2(:, i) - gc2
  end do

  !! set the H matrix
  matrix = matmul(coords1, transpose(coords2))

  !! call svd routine
  call svd(3, 3, matrix, u, smat, vt, ierr)
  if (ierr /= 0) then
     write (*, *) "error in svd for strucs:"
     write (*, *) nat1
     write (*, *) "struc1"
     do i = 1, nat1
        write (*, *) typ1(i), coords1(:, i)
     end do
     write (*, *) nat2
     write (*, *) "struc2"
     do i = 1, nat2
        write (*, *) typ2(i), coords2(:, i)
     end do
     write (*, *) "With the H matrix:"
     do i = 1, 3
        write (*, '(3f12.6)') matrix(i, :)
     end do
     flush (5)
     return
  end if

  !! set final rotation
  rmat = matmul(u, vt)

  !! set final translation
  translate = gc1 - matmul(rmat, gc2)

  !! test output
  ! write(*,*) 'translate'
  ! write(*,*) translate

  ! coords1(:,:) = coords1_in(:,:)
  ! coords2(:,:) = coords2_in(:,:)

  ! do i = 1, nat2
  !    coords2(:,i) = matmul(rmat,coords2(:,i)) + translate
  ! end do

  ! do i =1, nat2
  !    write(*,*) coords2(:,i)
  ! end do

end subroutine svdrot_m

!> @details
!! Optimal rotation by SVD.
!!
!! This routine forces the output rmat to be rotation (determinant = +1)
!!
!! The output can be applied as:
!!
!!~~~~~~~~~~~~~{.f90}
!!     coords2(:,i) = matmul( rmat, coords2(:,i) ) + translate
!!~~~~~~~~~~~~~
!!
!! OR
!!
!!~~~~~~~~~~~~~~~{.f90}
!!     coords1(:,i) = matmul( transpose(rmat), coords1(:,i) ) - &
!!                                  matmul( transpose(rmat), translate )
!!~~~~~~~~~~~~~~~
!!
!!
!! **Parameters:**
!! @param2{in, nat1, integer}    -> number of atoms in conf 1;
!! @param2{in, typ1, integer} size(nat1)    -> atomic types in conf 1;
!! @param2{in, coords1_in, real} size(3,nat1) -> coordinates of conf 1;
!! @param2{in, nat2, integer}    -> number of atoms in conf 2;
!! @param2{in, typ2, integer} size(nat2)    -> atomic types in conf 2;
!! @param2{in, coords2_in, real} size(3,nat2) -> coordinates of conf 2;
!! @param2{out, rmat, real} size(3,3)    -> 3x3 optimal rotation matrix
!! @param2{out, translate, real} size(3) -> optimal translation vector

!! @param[in] nat1    -> number of atoms in conf 1;
!! @param[in] typ1(nat1)    -> atomic types in conf 1;
!! @param[in] coords1_in(3,nat1) -> coordinates of conf 1;
!! @param[in] nat2    -> number of atoms in conf 2;
!! @param[in] typ2(nat2)    -> atomic types in conf 2;
!! @param[in] coords2_in(3,nat2) -> coordinates of conf 2;
!!
!! @param[out] rmat(3,3)    -> 3x3 optimal rotation matrix
!! @param[out] translate(3) -> optimal translation vector
!! @returns rmat, translate
!!
subroutine svd_forcerot(nat1, typ1, coords1_in, &
     nat2, typ2, coords2_in, &
     rmat, translate, ierr)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  real(rp), dimension(3, 3), intent(out) :: rmat
  real(rp), dimension(3), intent(out) :: translate
  integer(ip), intent(out) :: ierr

  real(rp), dimension(3, nat1) :: coords1
  real(rp), dimension(3, nat2) :: coords2
  real(rp), dimension(3) :: gc1, gc2
  real(rp), dimension(3, 3) :: matrix, u, smat, vt
  integer(ip) :: i
  real(rp) :: det_u, det_vt, det_s, det_m, det_r

  !! set initial copies
  coords1(:, :) = coords1_in(:, :)
  coords2(:, :) = coords2_in(:, :)

  !! set geo centers
  gc1(:) = 0.0_rp
  gc2(:) = 0.0_rp
  gc1(:) = sum(coords1(:, :), 2)/nat1
  gc2(:) = sum(coords2(:, :), 2)/nat2

  !! recenter
  do i = 1, nat1
     coords1(:, i) = coords1(:, i) - gc1
  end do
  do i = 1, nat2
     coords2(:, i) = coords2(:, i) - gc2
  end do

  !! set the H matrix
  matrix = matmul(coords1, transpose(coords2))
  ! write(*,*) 'Hmat'
  ! do i = 1, 3
  !    write(*,'(3(f12.8,x))') matrix(i,:)
  ! end do

  !! call svd routine
  call svd(3, 3, matrix, u, smat, vt, ierr)
  if (ierr /= 0) then
     write (*, *) "error in svd for strucs:"
     write (*, *) nat1
     write (*, *) "struc1"
     do i = 1, nat1
        write (*, *) typ1(i), coords1(:, i)
     end do
     write (*, *) nat2
     write (*, *) "struc2"
     do i = 1, nat2
        write (*, *) typ2(i), coords2(:, i)
     end do
     write (*, *) "With the H matrix:"
     do i = 1, 3
        write (*, '(3f12.6)') matrix(i, :)
     end do
     return
  end if

  ! call determinant3x3(u, det_u)
  ! write(*,*) 'u:', det_u
  ! do i = 1, 3
  !    write(*,*) u(i,:)
  ! end do

  ! call determinant3x3(smat, det_s)

  ! call determinant3x3(vt, det_vt)
  ! write(*,*) 'vt:',det_vt
  ! do i = 1, 3
  !    write(*,'(3f9.4)') vt(i,:)
  ! end do

  ! write(*,*) 'detm',det_m
  ! write(*,*) 'det_u',det_u
  ! write(*,*) 'det_s',det_s
  ! write(*,*) 'det_vt',det_vt

  !! set final rotation
  rmat = matmul(u, vt)
  call determinant3x3(rmat, det_r)

  !! force rmat to be rotation: Vt is written row-wise
  if (det_r .lt. -0.5_rp) vt(3, :) = -vt(3, :)
  rmat = matmul(u, vt)

  ! write(*,*) 'rmat:', det_r
  ! do i = 1, 3
  !    write(*,'(3f9.4)') rmat(i,:)
  ! end do

  !! set final translation
  translate = gc1 - matmul(rmat, gc2)

  !! test output
  ! write(*,*) 'translate'
  ! write(*,*) translate

  ! coords1(:,:) = coords1_in(:,:)
  ! coords2(:,:) = coords2_in(:,:)

  ! do i = 1, nat2
  !    coords2(:,i) = matmul(rmat,coords2(:,i)) + translate
  ! end do

  ! do i =1, nat2
  !    write(*,*) coords2(:,i)
  ! end do

end subroutine svd_forcerot

!> @details
!! Routine that does the loop over coords2 to find U_J, here called gamma.
!! This is the main IRA loop over the space of rotations.
!! Gamma basis can be returned as 3x3 matrix "gamma", or 3x1 integer
!! vector "gamma_idx", which contains indices of atoms that set up the
!! "gamma" basis. The first element gamma_idx(1) gives the reflection.
!!
!!
!! @param[in] nat1    -> number of atoms in conf 1;
!! @param[in] typ1_in(nat1)    -> atomic types in conf 1;
!! @param[in] coords1_in(3,nat1) -> coordinates of conf 1;
!! @param[in] nat2    -> number of atoms in conf 2;
!! @param[in] typ2_in(nat2)    -> atomic types in conf 2;
!! @param[in] coords2_in(3,nat2) -> coordinates of conf 2;
!! @param[in] kmax    -> distance cutoff for atoms included in search;
!!
!! @param[out] gamma(3,3)   -> The U_J reference frame (here called gamma), as 3x3 matrix;
!! @param[out] m_fin   -> Value giving reflection:
!!                            m_fin = 1 -> no reflection,
!!                            m_fin = -1 -> reflection.
!! @param[out] hd_out  -> Hausdorff distance value
!! @param[out] gamma_idx(3) -> possible intent(out) array of atomic indices giving gamma,
!!                                gamma_idx(1) = m_fin
!!
!! @param[inout] some_thr  -> threshold for dh, updated in self-sufficient way
!!
subroutine get_gamma_m(nat1, typ1_in, coords1_in, &
     nat2, typ2_in, coords2_in, &
     kmax, gamma, gamma_idx, hd_out, some_thr)
  use ira_precision
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1_in
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2_in
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  real(rp), intent(in) :: kmax
  real(rp), dimension(3, 3), intent(out) :: gamma
  integer(ip), dimension(3), intent(out) :: gamma_idx
  real(rp), intent(out) :: hd_out
  real(rp), intent(inout) :: some_thr

  integer(ip) :: i, j, idx1, idx2, m, k
  real(rp) :: hd, hd_old
  real(rp), dimension(2, nat2) :: d_o
  real(rp), dimension(nat2) :: dists
  integer(ip), dimension(nat2) :: found
  integer(ip), dimension(nat1) :: typ1
  real(rp), dimension(3, nat1) :: coords1
  integer(ip), dimension(nat2) :: typ2
  real(rp), dimension(3, nat2) :: coords2
  logical :: fail
  integer(ip) :: count
  integer(ip) :: m_fin

  !! set local copies
  typ1(:) = typ1_in(:)
  coords1(:, :) = coords1_in(:, :)
  typ2(:) = typ2_in(:)
  coords2(:, :) = coords2_in(:, :)
  !!
  !! sort coords2 by size
  !!
  do i = 1, nat2
     d_o(1, i) = norm2(coords2(:, i))
     d_o(2, i) = real(i, rp)
  end do
  !!
  call sort(nat2, 2, d_o, 1)
  !!
  !! permute to this order
  !!
  ! coords2(:,:) = coords2(:,nint(d_o(2,:)) )
  ! typ2(:) = typ2( nint(d_o(2,:)) )
  call permute_real_2d(nat2, 3, coords2, nint(d_o(2, :)))
  call permute_int_1d(nat2, typ2, nint(d_o(2, :)))
  !!
  !! search for gamma
  !!
  hd_old = 999.9_rp
  idx1 = 0
  idx2 = 0
  m_fin = 0
  !!
  count = 0
  do i = 1, nat2
     !!
     !! if out of range, exit this i
     !!
     if (norm2(coords2(:, i)) .gt. kmax) exit
     !!
     do j = 1, nat2
        !!
        if (i .eq. j) cycle
        !!
        !! if out of range, exit this j
        !!
        if (norm2(coords2(:, j)) .gt. kmax) exit
        !!
        !! set current bas
        !!
        m = 1
        call set_orthonorm_bas(coords2(:, i), coords2(:, j), gamma, fail)
        if (fail) cycle
        !!
        count = count + 1
        !!
        !! rotate coords_ev
        !!
        do k = 1, nat1
           coords1(:, k) = matmul(transpose(gamma), coords1(:, k))
        end do
        !!
        !! calc dists
        !!
        found(:) = 0_ip
        dists(:) = 0.0_rp
        call cshda(nat1, typ1, coords1, &
             nat2, typ2, coords2, some_thr, found, dists)
        !!
        !! Hausdorff
        hd = maxval(dists(1:nat1))
        !!
        !! keep the hd thr at minimum
        if (hd .lt. some_thr) some_thr = hd
        !!
        ! write(*,*) nat1+nat2
        ! write(*,*) 'bas',i,j, hd
        ! do k = 1, nat1
        !    write(*,*) 1, coords1(:,k)
        ! end do
        ! do k = 1, nat2
        !    write(*,*) 2, coords2(:,k)
        ! end do

        ! write(*,'(5i3,3f13.7)') i,j,nint(d_o(2,i)), nint(d_o(2,j)),m,hd, norm2(coords2(:,j))
        if (hd .lt. hd_old) then
           hd_old = hd
           idx1 = i
           idx2 = j
           m_fin = m
        end if
        !!
        !! rotate back
        !!
        do k = 1, nat1
           coords1(:, k) = matmul(gamma, coords1(:, k))
        end do
        !!
        !! do for m=-1 (reflection)
        !!
        m = -1
        gamma(3, :) = -gamma(3, :)
        !!
        !! rotate coords_ev
        !!
        do k = 1, nat1
           coords1(:, k) = matmul(transpose(gamma), coords1(:, k))
        end do
        !!
        !! calc dists
        !!
        call cshda(nat1, typ1, coords1, &
             nat2, typ2, coords2, some_thr, found, dists)
        !!
        !! Hausdorff
        hd = maxval(dists(1:nat1))
        !!
        !! keep the hd thr at minimum
        if (hd .lt. some_thr) some_thr = hd
        !!
        ! write(*,'(5i3,3f13.7)') i,j,nint(d_o(2,i)), nint(d_o(2,j)),m,hd, norm2(coords2(:,j))
        if (hd .lt. hd_old) then
           hd_old = hd
           idx1 = i
           idx2 = j
           m_fin = m
        end if

        !!
        !! rotate back
        !!
        do k = 1, nat1
           coords1(:, k) = matmul(gamma, coords1(:, k))
        end do
        !!
        !! early exit criterion idea
        ! if( hd_old .lt. some_threshold ) goto 111
     end do
  end do
  ! 111 continue
  !!
  ! write(*,'(a,x,f12.8)') 'hd', hd_old
  !!
  !! set found bas
  !!
  ! write(*,*) 'setting gama idx:',idx1,idx2,'m_fin:',m_fin, hd_old
  ! write(*,*) coords2(:,idx1)
  ! write(*,*) coords2(:,idx2)
  !!
  !! set data
  !!
  if (idx1 .eq. 0 .or. idx2 .eq. 0 .or. m_fin .eq. 0) then
     !! if nothing is found: BUG (probably too small kmax)
     !! But can also happen when searching nonequal nat...
     !! output gamma as identity matrix, all idx to zero
     gamma_idx(:) = 0
     gamma(:, :) = 0.0_rp
     gamma(1, 1) = 1.0_rp
     gamma(2, 2) = 1.0_rp
     gamma(3, 3) = 1.0_rp
  else
     !!
     !! set found result
     call set_orthonorm_bas(coords2(:, idx1), coords2(:, idx2), gamma, fail)
     gamma(3, :) = m_fin*gamma(3, :)

     !! for gamma_idx, use unsorted indices, which correspond to input order
     gamma_idx(1) = m_fin
     gamma_idx(2) = nint(d_o(2, idx1))
     gamma_idx(3) = nint(d_o(2, idx2))
  end if

  hd_out = hd_old

  ! write(*,*) 'nbas tested:', count
  ! write(*,*) 'gamma idx:',gamma_idx(:)
  return
end subroutine get_gamma_m

!> @details
!! Unified call to IRA for structures with equal and nonequal number of atoms.
!! The SVD is NOT performed at the end of this routine!
!!
!! the result is applied to struc 2:
!!
!! @code{.f90}
!!    j = permutation(i)
!!    coords2(:,i) = matmul( rotation, coords2(:,j) ) + tr
!! @endcode
!!
!! @param[in] nat1 :: number of atoms in struc 1
!! @param[in] typ1_in(nat1) :: atomic types of struc 1
!! @param[in] coords1_in(3,nat1) :: atoms positions of struc1
!! @param[in] candidate_1(nat1) :: candidate cenral atom of struc 1
!! @param[in] nat2 :: number of atoms of struc 2
!! @param[in] typ2_in(nat2) :: atomic types of struc 2
!! @param[in] coords2_in(3,nat2) :: atomic positions of struc 2
!! @param[in] candidate_2(nat2) :: candidate central atoms of struc 2
!! @param[in] kmax_factor :: multiplicative factor for kmax, should > 1.
!!
!! @param[out] rotation(3,3) :: R_apx rotation matrix after IRA
!! @param[out] translation(3) :: translation vector
!! @param[out] permutation(nat2) :: atomic permutations
!! @param[out] hd_out :: hausdorff distance
!! @param[out] ierr :: error value, negative on error, zero otherwise
subroutine ira_unify(nat1, typ1_in, coords1_in, candidate_1, &
     nat2, typ2_in, coords2_in, candidate_2, &
     kmax_factor, rotation, translation, permutation, hd_out, ierr)

  use ira_precision
  use err_module
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1_in
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), dimension(nat1), intent(in) :: candidate_1
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2_in
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  integer(ip), dimension(nat2), intent(in) :: candidate_2
  real(rp), intent(in) :: kmax_factor
  real(rp), dimension(3, 3), intent(out) :: rotation
  real(rp), dimension(3), intent(out) :: translation
  integer(ip), dimension(nat2), intent(out) :: permutation
  real(rp), intent(out) :: hd_out
  integer(ip), intent(out) :: ierr

  integer(ip) :: i, ii, j, jj
  integer(ip), dimension(nat1) :: typ1
  real(rp), dimension(3, nat1) :: coords1
  integer(ip), dimension(nat2) :: typ2
  real(rp), dimension(3, nat2) :: coords2
  integer(ip) :: c1, c2, c1min, c2min, idxm, idx1, idx2
  integer(ip), dimension(nat2) :: found
  real(rp), dimension(nat2) :: dists
  real(rp), dimension(2, nat1) :: d_o
  real(rp) :: hd_old, some_thr, hd
  real(rp), dimension(3) :: rc1, rc2
  real(rp), dimension(3, 3) :: beta, gamma, beta_min, gamma_min, invb
  logical :: fail
  real(rp) :: dist_k
  integer(ip), dimension(3) :: gamma_idx
  integer(ip) :: count

  ierr = ERR_OTHER

  !!
  !! REQUIREMENT: nat1 .le. nat2
  if (nat1 .gt. nat2) then
     write (*, *) "error in ira_unify: nat1 > nat2", nat1, nat2
     return
  end if

  !! make local copies of input structures
  typ1(:) = typ1_in(:)
  coords1(:, :) = coords1_in(:, :)
  typ2(:) = typ2_in(:)
  coords2(:, :) = coords2_in(:, :)

  !! initialize some variables
  c1min = 0
  c2min = 0
  idxm = 0
  idx1 = 0
  idx2 = 0
  beta_min(:, :) = 0.0_rp
  gamma_min(:, :) = 0.0_rp
  do i = 1, 3
     beta_min(i, i) = 1.0_rp
     gamma_min(i, i) = 1.0_rp
  end do
  permutation(:) = 0
  do i = 1, nat2
     permutation(i) = i
  end do

  rotation(:, :) = 0.0_rp
  rotation(1, 1) = 1.0_rp
  rotation(2, 2) = 1.0_rp
  rotation(3, 3) = 1.0_rp
  translation(:) = 0.0_rp
  hd_out = 999.9_rp

  hd_old = 999.8_rp
  some_thr = 9999.9_rp
  count = 0
  !!
  !! for each candidate center in struc 1:
  !!
  do ii = 1, nat1

     !! candidate index
     c1 = candidate_1(ii)

     !! we have exhausted all candidates in 1
     if (c1 .eq. 0) exit

     !!
     !! shift struc 1 to respective vector
     !!
     call select_rc(nat1, coords1, c1, rc1)
     !!
     do i = 1, nat1
        coords1(:, i) = coords1(:, i) - rc1
     end do

     !! sort by size (distance from center)
     do i = 1, nat1
        d_o(1, i) = norm2(coords1(:, i))
        d_o(2, i) = real(i)
     end do

     call sort(nat1, 2, d_o, 1)

     !! permute coords1 to that order
     ! coords1(:,:) = coords1(:,nint(d_o(2,:)))
     ! typ1(:) = typ1(nint(d_o(2,:)) )
     call permute_real_2d(nat1, 3, coords1, nint(d_o(2, :)))
     call permute_int_1d(nat1, typ1, nint(d_o(2, :)))

     !!
     !! Find some basis in structure 1. Originally we take the first possible
     !! basis that is found. If you have a better guess, modify this part.
     !! Ideas on optimal basis: the bases beta and gamma should be as equivalent
     !! as possible, therefore the atoms which set them should have as low
     !! distortions as possible.
     !!
     do i = 1, nat1
        do j = 1, nat1
           if (i .eq. j) cycle
           !! set bas
           ! write(*,*) coords1(:,i)
           ! write(*,*) coords1(:,j)
           ! write(*,*) i,j
           call set_orthonorm_bas(coords1(:, i), coords1(:, j), beta, fail)
           !! exit on first non-failed basis
           if (.not. fail) exit
        end do
        if (.not. fail) exit
     end do
     ! write(*,*) '1 in',i,j
     ! write(*,*) i, coords1(:,i)
     ! write(*,*) j, coords1(:,j)
     ! write(*,*)
     !!
     !! Measure "cutoff" of basis, multiply by a factor which is 1.2 by default.
     !! This factor should be greather than 1 in all cases. Decreasing it makes
     !! the search space for basis in structure 2 smaller, which can make the
     !! algorithm faster when structures are quite similar. It can also make it
     !! miss the optimal basis, in cases when the distortion is exactly on atoms
     !! corresponding to those which set the basis in structure 1.
     !!
     dist_k = max(norm2(coords1(:, i)), norm2(coords1(:, j)))
     dist_k = dist_k*kmax_factor
     !!
     !! rotate 1 to that bas
     !!
     do i = 1, nat1
        coords1(:, i) = matmul(beta, coords1(:, i))
     end do

     !!
     !! for each candidate center in struc 2:
     !!
     do jj = 1, nat2

        !! candidate index
        c2 = candidate_2(jj)

        !! we have exhausted all candidates in 2
        if (c2 .eq. 0) exit

        !!
        !! shift struc 2 to respective vector
        !!
        call select_rc(nat2, coords2, c2, rc2)
        !!
        do i = 1, nat2
           coords2(:, i) = coords2(:, i) - rc2
        end do

        !! find stuffs
        !!
        !! get gamma for this central atm
        !!
        call get_gamma_m(nat1, typ1(1:nat1), coords1(1:3, 1:nat1), &
             nat2, typ2(1:nat2), coords2(1:3, 1:nat2), &
             dist_k, gamma, gamma_idx, hd, some_thr)
        !!
        !! rotate to found basis
        !!
        do i = 1, nat1
           coords1(:, i) = matmul(transpose(gamma), coords1(:, i))
        end do

        !! if gamma is found:
        if (gamma_idx(1) .ne. 0) then
           !!
           count = count + 1
           !!
           !! find permutations
           !!
           call cshda(nat1, typ1(1:nat1), coords1(1:3, 1:nat1), &
                nat2, typ2(1:nat2), coords2(1:3, 1:nat2), &
                999.9_rp, found(1:nat2), dists(1:nat2))
        else
           !! all dists big
           dists(:) = 999.9_rp
        end if

        !! Hausdorff of this c1, up to nat1
        hd = maxval(dists(1:nat1))

        ! write(*,*) 'aa:',c1, c2, hd, some_thr, hd_old
        ! write(*,*) "gamma:"
        ! write(*, '(3f8.4)') gamma(1,:)
        ! write(*, '(3f8.4)') gamma(2,:)
        ! write(*, '(3f8.4)') gamma(3,:)
        !!
        if (hd .lt. hd_old) then
           hd_old = hd
           c1min = c1
           c2min = c2
           idxm = gamma_idx(1)
           idx1 = gamma_idx(2)
           idx2 = gamma_idx(3)
           beta_min = beta
           gamma_min = gamma
        end if

        !! rotate struc 1 back to orig
        do i = 1, nat1
           coords1(:, i) = matmul(gamma, coords1(:, i))
        end do

        !! shift struc 2 back
        do i = 1, nat2
           coords2(:, i) = coords2(:, i) + rc2
        end do

     end do

     !!
     !! rotate 1 back to original
     !!
     do i = 1, nat1
        coords1(:, i) = matmul(transpose(beta), coords1(:, i))
     end do

     !! shift struc 1 back
     do i = 1, nat1
        coords1(:, i) = coords1(:, i) + rc1
     end do

     !! permute struc1 back to orig
     ! coords1(:,nint(d_o(2,:))) = coords1(:,:)
     ! typ1(nint(d_o(2,:)) ) = typ1(:)
     call permute_real_2d_back(nat1, 3, coords1, nint(d_o(2, :)))
     call permute_int_1d_back(nat1, typ1, nint(d_o(2, :)))

  end do

  ! write(*,*) "beta_min:"
  ! write(*,'(3f8.4)') beta_min(1,:)
  ! write(*,'(3f8.4)') beta_min(2,:)
  ! write(*,'(3f8.4)') beta_min(3,:)
  ! write(*,*) "gamma_min:"
  ! write(*,'(3f8.4)') gamma_min(1,:)
  ! write(*,'(3f8.4)') gamma_min(2,:)
  ! write(*,'(3f8.4)') gamma_min(3,:)

  !!
  !! at the end of this loop, the structures 1 and 2 should be
  !! identical to the input structures
  !!

  !! no attempts have been made, return error (too small dist_k?)
  if (count .eq. 0) then
     ierr = ERR_TOO_SMALL_KMAX
     ! write(*,*) "count=0", ierr
     return
  end if

  !! nothing has been found, this is probably same reason: too small dist_k
  if (c1min .eq. 0 .or. c2min .eq. 0) then
     ierr = ERR_OTHER
     ! write(*,*) "another err"
     return
  end if

  !!
  !! set found data
  !!
  !! center vectors
  call select_rc(nat1, coords1, c1min, rc1)
  call select_rc(nat2, coords2, c2min, rc2)
  !!
  !! rotation
  invb = transpose(beta_min)
  rotation = matmul(invb, gamma_min)
  !!
  !! translation
  translation(:) = rc1 - matmul(rotation, rc2)

  !! rotate + translate struc2
  do i = 1, nat2
     coords2(:, i) = matmul(rotation, coords2(:, i)) + translation
  end do

  !! permutation
  if (idxm .ne. 0) then
     !! if no bug in gamma (which happens if dist_k too small):
     !! find final permutations
     call cshda(nat1, typ1, coords1, &
          nat2, typ2, coords2, &
          999.9_rp, found, dists)
     !!
     !! set final permutation
     !!
     permutation(:) = found(:)
  else
     !! all dists large
     dists(:) = 999.9_rp
  end if

  !! set hd
  hd_out = maxval(dists(1:nat1))
  
  ierr = 0

end subroutine ira_unify

!> @details
!! set_candidates + ira_unify + svd
!! candidates are set as per default:
!! <ul class="no-bullets">
!!   <li> gc if equal; </li>
!!   <li> 1st atom in smaller, </li>
!!   <li> all in larger </li>
!!   <li> <code> hah</code>
!! </ul>
!!
!! **Input**:
!! @param2{in, nat1, integer    }               :: number of atoms in struc 1
!! @param2{in, typ1_in, integer } size(nat1)    :: atomic types of struc 1
!! @param2{in, coords1_in, real } size(3, nat1) :: atomic positions of struc 1
!! @param2{in, nat2, integer    }               :: number of atoms of struc 2
!! @param2{in, typ2_in, integer } size(nat2)    :: atomic types of struc 2
!! @param2{in, coords2_in, real } size(3, nat2) :: atomic positions of struc 2
!! @param2{in, kmax_factor, real }              :: multiplicative factor for kmax, should > 1.
!!
!! **Output**:
!! @param2{out, rotation, real } size(3,3)      :: \f$ R_{apx}\f$ rotation matrix after IRA
!! @param2{out, translation, real } size(3)     :: translation vector
!! @param2{out, permutation, integer } size(nat2) :: atomic permutations
!! @param2{out, hd, real}                       :: hausdorff distance
!! @param2{out, rmsd, real }                    :: rmsd distance
!! @param2{out, ierr, integer }                 :: error value, negative on error, zero otherwise
!!
!!

!! @param[in] <code>nat1</code> :: number of atoms in struc 1
!! @param[in] `typ1_in(nat1)` :: atomic types of struc 1
!! @param[in] coords1_in(3,nat1) :: atoms positions of struc1
!! @param[in] nat2 :: number of atoms of struc 2
!! @param[in] typ2_in(nat2) :: atomic types of struc 2
!! @param[in] coords2_in(3,nat2) :: atomic positions of struc 2
!! @param[in] kmax_factor :: multiplicative factor for kmax, should > 1.
!!
!! @param[out] rotation(3,3) :: R_apx rotation matrix after IRA
!! @param[out] translation(3) :: translation vector
!! @param[out] permutation(nat2) :: atomic permutations
!! @param[out] hd :: hausdorff distance
!! @param[out] rmsd :: rmsd distance
!! @param[out] ierr :: error value, negative on error, zero otherwise
subroutine ira_svd(nat1, typ1_in, coords1_in, &
     nat2, typ2_in, coords2_in, &
     kmax_factor, rotation, translation, permutation, &
     hd, rmsd, ierr)
  use ira_precision
  use err_module, only: get_err_msg
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1_in
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2_in
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  real(rp), intent(in) :: kmax_factor
  real(rp), dimension(3, 3), intent(out) :: rotation
  real(rp), dimension(3), intent(out) :: translation
  integer(ip), dimension(nat2), intent(out) :: permutation
  real(rp), intent(out) :: hd
  real(rp), intent(out) :: rmsd
  integer(ip), intent(out) :: ierr

  integer(ip), dimension(nat1) :: typ1
  real(rp), dimension(3, nat1) :: coords1
  integer(ip), dimension(nat2) :: typ2
  real(rp), dimension(3, nat2) :: coords2
  integer(ip), dimension(nat1) :: candidate_1
  integer(ip), dimension(nat2) :: candidate_2
  integer(ip) :: i
  real(rp) :: hd_out
  real(rp), dimension(3, 3) :: svd_rot
  real(rp), dimension(3) :: rdum, svd_tr

  !! copy input data
  typ1(:) = typ1_in(:)
  typ2(:) = typ2_in(:)
  coords1(:, :) = coords1_in(:, :)
  coords2(:, :) = coords2_in(:, :)

  !!
  !! form candidates for central atms in 1 and 2
  !!
  call set_candidates(nat1, typ1, coords1, &
       nat2, typ2, coords2, &
       candidate_1, candidate_2)

  !!
  !! call main ira with candidates
  !!
  call ira_unify(nat1, typ1, coords1, candidate_1, &
       nat2, typ2, coords2, candidate_2, &
       kmax_factor, rotation, translation, permutation, hd_out, ierr)
  if (ierr .ne. 0) then
     ! write(*,*) "error in ira_unify! ierr code:", ierr
     ! write(*,*) get_err_msg( ierr )
     return
  end if
  !!
  !! apply found transformation
  !!
  typ2(:) = typ2(permutation(:))
  coords2(:, :) = coords2(:, permutation(:))
  !!
  do i = 1, nat2
     coords2(:, i) = matmul(rotation, coords2(:, i)) + translation
  end do

  !!
  !! call SVD
  !!
  call svdrot_m(nat1, typ1, coords1, &
       nat1, typ2(1:nat1), coords2(:, 1:nat1), &
       svd_rot, svd_tr, ierr)
  if (ierr /= 0) then
     return
  end if

  !!
  !! apply svd
  !!
  ! do i = 1, nat2
  !    coords2(:,i) = matmul( svd_rot, coords2(:,i) ) + svd_tr
  ! end do

  !! rotate back to orig
  do i = 1, nat2
     coords2(:, i) = matmul(transpose(rotation), coords2(:, i)) - matmul(transpose(rotation), translation)
  end do

  !!
  !! put together apx and svd
  !!
  rotation = matmul(svd_rot, rotation)
  translation = matmul(svd_rot, translation) + svd_tr
  !! apply
  do i = 1, nat2
     coords2(:, i) = matmul(rotation, coords2(:, i)) + translation
  end do
  ! do i = 1, nat1
  !    coords1(:,i) = matmul(transpose(rotation),coords1(:,i)) - matmul(transpose(rotation),translation)
  ! end do

  !!
  !! compute dH and rmsd distances for output
  !!
  hd = 0.0_rp
  do i = 1, nat1
     hd = max(hd, norm2(coords1(:, i) - coords2(:, i)))
  end do

  rmsd = 0.0_rp
  do i = 1, nat1
     rdum = coords1(:, i) - coords2(:, i)
     rmsd = rmsd + dot_product(rdum, rdum)
  end do
  rmsd = sqrt(rmsd/nat1)

end subroutine ira_svd

subroutine ira_get_err_msg(ierr, msg)
  use ira_precision
  use err_module, only: get_err_msg
  implicit none
  integer(ip), intent(in) :: ierr
  character(512), intent(out) :: msg

  character(:), allocatable :: me

  allocate (me, source=get_err_msg(ierr))
  msg = me
  deallocate (me)
end subroutine ira_get_err_msg
