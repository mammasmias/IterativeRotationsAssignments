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


  ! subroutine determinant( a, d )
  !   !! determinant of a 3x3 matrix a
  !   implicit none
  !   real, dimension(3,3), intent(in) :: a
  !   real, intent(out) :: d

  !   d = a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) + &
  !         a(1,2)*a(2,3)*a(3,1) - a(2,2)*a(3,1)*a(1,3) + &
  !         a(2,1)*a(3,2)*a(1,3) - a(3,3)*a(1,2)*a(2,1)

  ! end subroutine determinant


  subroutine sort( n, ndim, array, axis )
    !! Use mergesort to sort the input array by the chosen axis.
    !! The routine mergesort is located in sorting_module.
    use sorting_module, only: mergesort
    implicit none
    integer,                 intent(in) :: n
    integer,                 intent(in) :: ndim
    real, dimension(ndim,n), intent(inout) :: array
    integer,                 intent(in) :: axis

    real, allocatable :: work(:,:)

    !! allocate working space for the sorting routine
    allocate( work(1:ndim, 1:(n+1)/2), source = 0.0)

    call mergesort( array, work, axis )

    deallocate( work )

  end subroutine sort


  subroutine svd( m, n, a, u, s, v )
    !! routine that calls LAPACK svd routine
    !!=========================================
    !!
    !! intent(in):
    !! m     -> leading dimension of matrix a
    !! n     -> second dimension of matrix a
    !! a     -> matrix a
    !!
    !! intent(out):
    !! s     -> diagonal matrix of singular values
    !! u     -> orthonormal matrix u
    !! v     -> orthonormal matrix v
    implicit none

    integer,              intent(in) :: m
    integer,              intent(in) :: n
    real, dimension(m,n), intent(in) :: a
    real, dimension(m,n), intent(out) :: s
    real, dimension(m,n), intent(out) :: u
    real, dimension(m,n), intent(out) :: v

    real, dimension(m,n) :: a_copy
    real, dimension(min(m,n)) :: sdiag
    real, allocatable :: work(:)
    integer :: lwork, i, info, lda, ldu, ldv
    character(len=3) :: jobu,jobv


    lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

    allocate ( work(1:lwork) )
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
    a_copy(1:m,1:n) = a(1:m,1:n)

    !!
    !! for the single precision real use this:
    !call sgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
    !     lwork, info )
    !!
    !!
    !! for the double precision real use this:
    call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
         lwork, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The SVD could not be calculated.'
      write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
      write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
      return
    end if
    !
    !  Make the MxN matrix S from the diagonal values in SDIAG.
    !
    s(1:m,1:n) = 0.0D+00
    do i = 1, min ( m, n )
      s(i,i) = sdiag(i)
    end do
    !
    !  Transpose V.
    !
    !    v = transpose ( v )

    deallocate ( work )

    return

  end subroutine svd


  ! subroutine cross_prod( a, b, c )
  !   !> @brief Cross product of two vectors
  !   implicit none
  !   real, dimension(3), intent(in) :: a
  !   real, dimension(3), intent(in) :: b
  !   real, dimension(3), intent(out) :: c

  !   c(1) = a(2)*b(3) - a(3)*b(2)
  !   c(2) = a(3)*b(1) - a(1)*b(3)
  !   c(3) = a(1)*b(2) - a(2)*b(1)

  ! end subroutine cross_prod



  subroutine permute_real_2d( n, m, array, order )
    !! permute a real 2D array into order, equivalent to:
    !!
    !!   array(:,:) = array(:, order(:) )
    !!
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    real, dimension(m, n), intent(inout) :: array
    integer, dimension(n), intent(in) :: order

    integer :: i
    real, allocatable :: tmp(:,:)

    !! tmp copy
    allocate( tmp(1:m, 1:n), source=array )

    !! permute
    do i = 1, n
       array(:,i) = tmp(:, order(i) )
    end do

    deallocate( tmp )

  end subroutine permute_real_2d

  subroutine permute_real_2d_back( n, m, array, order )
    !! permute a real 2D array into inverse order,
    !! equivalent to:
    !!
    !!   array(:, order(:) ) = array(:,:)
    !!
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    real, dimension(m, n), intent(inout) :: array
    integer, dimension(n), intent(in) :: order

    integer :: i
    real, allocatable :: tmp(:,:)

    !! tmp copy
    allocate( tmp(1:m, 1:n), source=array )

    !! permute
    do i = 1, n
       array(:, order(i)) = tmp(:, i )
    end do

    deallocate( tmp )

  end subroutine permute_real_2d_back


  subroutine permute_int_1d( n, array, order )
    !! permute an integer 1D array into order, equivalent to:
    !!
    !!    array(:) = array( order(:) )
    !!
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(inout) :: array
    integer, dimension(n), intent(in) :: order

    integer :: i
    integer, allocatable :: tmp(:)

    !! tmp copy
    allocate( tmp(1:n), source=array )

    !! permute
    do i = 1, n
       array(i) = tmp( order(i) )
    end do

    deallocate( tmp )

  end subroutine permute_int_1d

  subroutine permute_int_1d_back( n, array, order )
    !! permute an integer 1D array into inverse order,
    !! equivalent to:
    !!
    !!    array( order(:) ) = array(:)
    !!
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(inout) :: array
    integer, dimension(n), intent(in) :: order

    integer :: i
    integer, allocatable :: tmp(:)

    !! tmp copy
    allocate( tmp(1:n), source=array )

    !! permute
    do i = 1, n
       array( order(i) ) = tmp( i )
    end do

    deallocate( tmp )

  end subroutine permute_int_1d_back


  subroutine set_orthonorm_bas( vec1, vec2, basis, fail )
    !> @brief Set orthonormal basis from input vectors, third is cross of
    !! the first two, fail happens when input is collinear.
    !!
    !! basis on output contains basis vectors in rows.
    implicit none

    real, dimension(3),   intent(in) :: vec1, vec2
    real, dimension(3,3), intent(out) :: basis
    logical,              intent(out) :: fail

    real :: prod
    real :: collinearity_thr
    real :: small_size_thr
    real :: norm_v

    !! threshold for collinearity of two vectors
    collinearity_thr = 0.95

    !! threshold for too small vectors
    !! (for numerical reasons, don't want to normalize a too small vector)
    small_size_thr = 1e-2

    fail = .true.
    basis(:,:) = 0.0
    ! write(*,*) vec1
    ! write(*,*) vec2

    !! first vector, normalize
    norm_v = sqrt( dot_product(vec1,vec1) )
    basis(1,:) = vec1(:) / norm_v
    !! if too small, return with fail = .true.
    if( norm_v .lt. small_size_thr ) return

    !! second vector, normalize
    norm_v = sqrt( dot_product(vec2,vec2) )
    basis(2,:) = vec2(:) / norm_v
    !! if too small, return with fail = .true.
    if( norm_v .lt. small_size_thr ) return

    !! check projection
    prod = dot_product( basis(1,:), basis(2,:) )

    !! if vectors are collinear, return with fail = .true.
    if( abs(prod) .gt. collinearity_thr ) then
       fail = .true.
       ! write(*,*) 'failed collinearity',abs(prod)
       return
    endif

    !! take orthogonal component of second vector, normalize
    basis(2,:) = basis(2,:) - prod*basis(1,:)
    norm_v = sqrt( dot_product(basis(2,:), basis(2,:)) )
    basis(2,:) = basis(2,:) / norm_v

    !! third vector is cross product of first two, normalize
    basis(3,1) = basis(1,2)*basis(2,3) - basis(1,3)*basis(2,2)
    basis(3,2) = basis(1,3)*basis(2,1) - basis(1,1)*basis(2,3)
    basis(3,3) = basis(1,1)*basis(2,2) - basis(1,2)*basis(2,1)

    norm_v = sqrt( dot_product(basis(3,:), basis(3,:)) )
    basis(3,:) = basis(3,:) / norm_v

    !!
    !! check for NaN. Can it happen?
    !!
    if( any(basis .ne. basis) ) then
       fail = .true.
       ! write(*,*) 'fail NaN'
       ! write(*,*) basis(1,:)
       ! write(*,*) basis(2,:)
       ! write(*,*) basis(3,:)
       return
    endif

    !! if we come to this point, the basis is ok
    fail = .false.
    return
  end subroutine set_orthonorm_bas


  subroutine svdrot_m( nat1, typ1, coords1_in, &
                       nat2, typ2, coords2_in, &
                       rmat, translate)
    !> @detail
    !! Optimal rotation by SVD.
    !!
    !! This rouine allows rotation matrix to be a reflection.
    !! no 'correction' when determinant is negative
    !!
    !! The output can be applied as:
    !!
    !!     coords2(:,i) = matmul( rmat, coords2(:,i) ) + translate
    !!
    !! OR
    !!
    !!     coords1(:,i) = matmul( transpose(rmat), coords1(:,i) ) - &
    !!                                  matmul( transpose(rmat), translate )
    !!
    !!================================================
    !!
    !! intent(in):
    !! nat1    -> number of atoms in conf 1;
    !! typ1    -> atomic types in conf 1;
    !! coords1 -> coordinates of conf 1;
    !! nat2    -> number of atoms in conf 2;
    !! typ2    -> atomic types in conf 2;
    !! coords2 -> coordinates of conf 2;
    !!
    !! intent(out):
    !! rmat    -> 3x3 optimal rotation matrix
    !! translate -> optimal translation vector
    !!
    implicit none
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1
    real, dimension(3,nat1),  intent(in) :: coords1_in
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2
    real, dimension(3,nat2),  intent(in) :: coords2_in
    real,dimension(3,3),      intent(out) :: rmat
    real, dimension(3),       intent(out) :: translate

    real, allocatable :: coords1(:,:)
    real, allocatable :: coords2(:,:)
    real, dimension(3) :: gc1, gc2
    real, dimension(3,3) :: matrix, u, smat, vt
    integer :: i

    !! set initial copies
    allocate( coords1(1:3,1:nat1), source=coords1_in )
    allocate( coords2(1:3,1:nat2), source=coords2_in )

    !! set geo centers
    gc1(:) = 0.0
    gc2(:) = 0.0
    gc1(:) = sum( coords1(:,:), 2 )/nat1
    gc2(:) = sum( coords2(:,:), 2 )/nat2

    !! recenter
    do i = 1, nat1
       coords1(:,i) = coords1(:,i) - gc1
    end do
    do i = 1, nat2
       coords2(:,i) = coords2(:,i) - gc2
    end do

    !! set the H matrix
    matrix = matmul( coords1, transpose(coords2))

    !! call svd routine
    call svd(3, 3, matrix, u, smat, vt )

    !! set final rotation
    rmat = matmul(u, vt)

    !! set final translation
    translate = gc1 - matmul(rmat,gc2)

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

    deallocate( coords1, coords2 )

  end subroutine svdrot_m


  subroutine get_gamma_m(nat1, typ1_in, coords1_in, &
                         nat2, typ2_in, coords2_in, &
                         kmax, gamma, gamma_idx, hd_out, some_thr )
    !> @detail
    !! Routine that does the loop over coords2 to find U_J, here called gamma.
    !! This is the main IRA loop over the space of rotations.
    !! Gamma basis can be returned as 3x3 matrix "gamma", or 3x1 integer
    !! vector "gamma_idx", which contains indices of atoms that set up the
    !! "gamma" basis. The first element gamma_idx(1) gives the reflection.
    !!
    !!================================================
    !! intent(in):
    !! nat1    -> number of atoms in conf 1;
    !! typ1    -> atomic types in conf 1;
    !! coords1 -> coordinates of conf 1;
    !! nat2    -> number of atoms in conf 2;
    !! typ2    -> atomic types in conf 2;
    !! coords2 -> coordinates of conf 2;
    !! kmax    -> distance cutoff for atoms included in search;
    !!
    !! intent(out):
    !! gamma   -> The U_J reference frame (here called gamma), as 3x matrix;
    !! m_fin   -> Value giving reflection:
    !!                        m_fin = 1 -> no reflection,
    !!                        m_fin = -1 -> reflection.
    !! hd_out  -> Hausdorff distance value
    !! gamma_idx -> possible intent(out) array of atomic indices giving gamma,
    !!              gamma_idx(1) = m_fin
    !!
    !! intent(inout):
    !! some_thr  -> threshold for dh, updated in self-sufficient way
    !!
    implicit none
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1_in
    real, dimension(3,nat1),  intent(in) :: coords1_in
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2_in
    real, dimension(3,nat2),  intent(in) :: coords2_in
    real,                     intent(in) :: kmax
    real, dimension(3,3),     intent(out) :: gamma
    integer, dimension(3),    intent(out) :: gamma_idx
    real,                     intent(out) :: hd_out
    real,                     intent(inout) :: some_thr

    integer :: i, j, idx1, idx2, m, k
    real :: hd, hd_old
    real, allocatable :: dists(:), d_o(:,:)
    integer, allocatable :: found(:)
    integer, dimension(nat1) :: typ1
    real, dimension(3,nat1) :: coords1
    integer, dimension(nat2) :: typ2
    real, dimension(3,nat2) :: coords2
    logical :: fail
    integer :: count
    integer :: m_fin


    !! set local copies
    typ1(:) = typ1_in(:)
    coords1(:,:) = coords1_in(:,:)
    typ2(:) = typ2_in(:)
    coords2(:,:) = coords2_in(:,:)
    !!
    !! sort coords2 by size
    !!
    allocate( d_o(1:2, 1:nat2) )
    do i = 1, nat2
       d_o(1,i) = norm2( coords2(:,i))
       d_o(2,i) = real(i)
    end do
    !!
    call sort( nat2, 2, d_o, 1)
    !!
    !! permute to this order
    !!
    ! coords2(:,:) = coords2(:,nint(d_o(2,:)) )
    ! typ2(:) = typ2( nint(d_o(2,:)) )
    call permute_real_2d( nat2, 3, coords2, nint(d_o(2,:)) )
    call permute_int_1d( nat2, typ2, nint(d_o(2,:)) )
    !!
    !! search for gamma
    !!
    hd_old = 999.9
    allocate( found(1:nat2) )
    allocate( dists(1:nat2) )
    idx1 = 0
    idx2 = 0
    m_fin = 0
    !!
    count = 0
    do i = 1, nat2
       !!
       !! if out of range, exit this i
       !!
       if( norm2( coords2(:,i)) .gt. kmax ) exit
       !!
       do j = 1, nat2
          !!
          if( i .eq. j ) cycle
          !!
          !! if out of range, exit this j
          !!
          if( norm2(coords2(:,j)) .gt. kmax ) exit
          !!
          !! set current bas
          !!
          m = 1
          call set_orthonorm_bas( coords2(:,i), coords2(:,j), gamma, fail )
          if( fail ) cycle
          !!
          count = count + 1
          !!
          !! rotate coords_ev
          !!
          do k = 1, nat1
             coords1(:,k) = matmul( transpose(gamma), coords1(:,k) )
          end do
          !!
          !! calc dists
          !!
          found(:) = 0
          dists(:) = 0.0
          call cshda( nat1, typ1, coords1, &
               nat2, typ2, coords2, some_thr, found, dists )
          !!
          !! Hausdorff
          hd = maxval(dists(1:nat1) )
          !!
          !! keep the hd thr at minimum
          if( hd .lt. some_thr ) some_thr = hd
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
          if( hd .lt. hd_old ) then
             hd_old = hd
             idx1 = i
             idx2 = j
             m_fin = m
          endif
          !!
          !! rotate back
          !!
          do k = 1, nat1
             coords1(:,k) = matmul( gamma, coords1(:,k) )
          end do
          !!
          !! do for m=-1 (reflection)
          !!
          m = -1
          gamma(3,:) = -gamma(3,:)
          !!
          !! rotate coords_ev
          !!
          do k = 1, nat1
             coords1(:,k) = matmul( transpose(gamma), coords1(:,k) )
          end do
          !!
          !! calc dists
          !!
          call cshda( nat1, typ1, coords1, &
               nat2, typ2, coords2, some_thr, found, dists )
          !!
          !! Hausdorff
          hd = maxval(dists(1:nat1))
          !!
          !! keep the hd thr at minimum
          if( hd .lt. some_thr ) some_thr = hd
          !!
          ! write(*,'(5i3,3f13.7)') i,j,nint(d_o(2,i)), nint(d_o(2,j)),m,hd, norm2(coords2(:,j))
          if( hd .lt. hd_old ) then
             hd_old = hd
             idx1 = i
             idx2 = j
             m_fin = m
          endif

          !!
          !! rotate back
          !!
          do k = 1, nat1
             coords1(:,k) = matmul( gamma, coords1(:,k) )
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
    if( idx1 .eq. 0 .or. idx2 .eq. 0 .or. m_fin .eq. 0) then
       !! if nothing is found: BUG (probably too small kmax)
       !! output gamma as identity matrix, all idx to zero
       gamma_idx(:) = 0
       gamma(:,:) = 0.0
       gamma(1,1) = 1.0
       gamma(2,2) = 1.0
       gamma(3,3) = 1.0
    else
       !!
       !! set found result
       call set_orthonorm_bas( coords2(:,idx1), coords2(:,idx2), gamma, fail)
       gamma(3,:) = m_fin*gamma(3,:)

       !! for gamma_idx, use unsorted indices, which correspond to input order
       gamma_idx(1) = m_fin
       gamma_idx(2) = nint(d_o(2,idx1))
       gamma_idx(3) = nint(d_o(2,idx2))
    endif

    hd_out = hd_old
    deallocate( d_o )
    deallocate( found, dists )

    ! write(*,*) 'nbas tested:', count
    ! write(*,*) 'gamma idx:',gamma_idx(:)
    return
  end subroutine get_gamma_m


  subroutine ira_unify( nat1, typ1_in, coords1_in, candidate_1, &
                        nat2, typ2_in, coords2_in, candidate_2, &
                        kmax_factor, rotation, translation, permutation, hd_out )
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1_in
    real, dimension(3, nat1), intent(in) :: coords1_in
    integer, dimension(nat1), intent(in) :: candidate_1
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2_in
    real, dimension(3,nat2),  intent(in) :: coords2_in
    integer, dimension(nat2), intent(in) :: candidate_2
    real,                     intent(in) :: kmax_factor
    real, dimension(3,3),     intent(out) :: rotation
    real, dimension(3),       intent(out) :: translation
    integer, dimension(nat2), intent(out) :: permutation
    real,                     intent(out) :: hd_out

    integer :: i, ii, j, jj
    integer, dimension(nat1) :: typ1
    real, dimension(3,nat1) :: coords1
    integer, dimension(nat2) :: typ2
    real, dimension(3,nat2) :: coords2
    integer :: c1, c2, c1min, c2min, idxm, idx1, idx2
    integer, dimension(nat2) :: found
    real, dimension(nat2) :: dists
    real, dimension(2,nat1) :: d_o
    real :: hd_old, some_thr, hd
    real, dimension(3) :: rc1, rc2
    real, dimension(3,3) :: beta, gamma, beta_min, gamma_min, invb
    logical :: fail
    real :: dist_k
    integer, dimension(3) :: gamma_idx


    !!
    !! REQUIREMENT: nat1 .le. nat2
    if( nat1 .gt. nat2) then
      write(*,*) 'error in ira_unify: nat1 > nat2',nat1, nat2
      stop
    endif

    !! make local copies of input structures
    typ1(:) = typ1_in(:)
    coords1(:,:) = coords1_in(:,:)
    typ2(:) = typ2_in(:)
    coords2(:,:) = coords2_in(:,:)


    !! initialize some variables
    c1min = 0
    c2min = 0
    idxm = 0
    idx1 = 0
    idx2 = 0
    beta_min(:,:) = 0.0
    gamma_min(:,:) = 0.0
    do i = 1, 3
      beta_min(i,i) = 1.0
      gamma_min(i,i) = 1.0
    end do
    permutation(:) = 0
    do i = 1, nat2
      permutation(i) = i
    end do


    hd_old = 999.9
    some_thr = 9999.9
    !!
    !! for each candidate center in struc 1:
    !!
    do ii = 1, nat1

      !! candidate index
      c1 = candidate_1(ii)

      !! we have exhausted all candidates in 1
      if( c1 .eq. 0 ) exit

      !!
      !! shift struc 1 to respective vector
      !!
      call select_rc( nat1, coords1, c1, rc1 )
      !!
      do i = 1, nat1
          coords1(:,i) = coords1(:,i) - rc1
      end do

      !! sort by size (distance from center)
      do i = 1, nat1
          d_o(1,i) = norm2(coords1(:,i))
          d_o(2,i) = real(i)
      end do

      call sort(nat1, 2, d_o, 1)

      !! permute coords1 to that order
      ! coords1(:,:) = coords1(:,nint(d_o(2,:)))
      ! typ1(:) = typ1(nint(d_o(2,:)) )
      call permute_real_2d( nat1, 3, coords1, nint(d_o(2,:)) )
      call permute_int_1d( nat1, typ1, nint(d_o(2,:)) )



      !!
      !! Find some basis in structure 1. Originally we take the first possible
      !! basis that is found. If you have a better guess, modify this part.
      !! Ideas on optimal basis: the bases beta and gamma should be as equivalent
      !! as possible, therefore the atoms which set them should have as low
      !! distortions as possible.
      !!
      do i = 1, nat1
          do j = 1, nat1
            if( i .eq. j ) cycle
            !! set bas
            ! write(*,*) coords1(:,i)
            ! write(*,*) coords1(:,j)
            ! write(*,*) i,j
            call set_orthonorm_bas( coords1(:,i),coords1(:,j), beta, fail)
            !! exit on first non-failed basis
            if( .not. fail ) exit
          end do
          if( .not. fail ) exit
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
      dist_k = max( norm2(coords1(:,i)), norm2(coords1(:,j)) )
      dist_k = dist_k*kmax_factor
      !!
      !! rotate 1 to that bas
      !!
      do i = 1, nat1
          coords1(:,i) = matmul(beta, coords1(:,i))
      end do


      !!
      !! for each candidate center in struc 2:
      !!
      do jj = 1, nat2

          !! candidate index
          c2 = candidate_2(jj)

          !! we have exhausted all candidates in 2
          if( c2 .eq. 0 ) exit

          !!
          !! shift struc 2 to respective vector
          !!
          call select_rc( nat2, coords2, c2, rc2 )
          !!
          do i = 1, nat2
            coords2(:,i) = coords2(:,i) - rc2
          end do


          !! find stuffs
          !!
          !! get gamma for this central atm
          !!
          call get_gamma_m(nat1, typ1(1:nat1), coords1(1:3,1:nat1), &
              nat2, typ2(1:nat2), coords2(1:3,1:nat2), &
              dist_k, gamma, gamma_idx, hd, some_thr )


          !!
          !! rotate to found basis
          !!
          do i = 1, nat1
            coords1(:,i) = matmul(transpose(gamma),coords1(:,i))
          end do

          !! if gamma is found:
          if( gamma_idx(1) .ne. 0 ) then
            !!
            !! find permutations
            !!
            call cshda( nat1, typ1(1:nat1), coords1(1:3,1:nat1), &
                  nat2, typ2(1:nat2), coords2(1:3,1:nat2), &
                  999.9, found(1:nat2), dists(1:nat2) )
          else
            !! all dists big
            dists(:) = 999.9
          endif

          !! Hausdorff of this c1, up to nat1
          hd = maxval(dists(1:nat1))


          ! write(*,*) 'aa:',c1, c2, hd, some_thr, hd_old
          !!
          if( hd .lt. hd_old ) then
            hd_old = hd
            c1min = c1
            c2min = c2
            idxm = gamma_idx(1)
            idx1 = gamma_idx(2)
            idx2 = gamma_idx(3)
            beta_min = beta
            gamma_min = gamma
          endif

          !! rotate struc 1 back to orig
          do i = 1, nat1
            coords1(:,i) = matmul(gamma, coords1(:,i))
          end do

          !! shift struc 2 back
          do i = 1, nat2
            coords2(:,i) = coords2(:,i) + rc2
          end do

      end do

      !!
      !! rotate 1 back to original
      !!
      do i = 1, nat1
          coords1(:,i) = matmul(transpose(beta), coords1(:,i))
      end do

      !! shift struc 1 back
      do i = 1, nat1
          coords1(:,i) = coords1(:,i) + rc1
      end do

      !! permute struc1 back to orig
      ! coords1(:,nint(d_o(2,:))) = coords1(:,:)
      ! typ1(nint(d_o(2,:)) ) = typ1(:)
      call permute_real_2d_back( nat1, 3, coords1, nint(d_o(2,:)) )
      call permute_int_1d_back( nat1, typ1, nint(d_o(2,:)) )


    end do

    !!
    !! at the end of this loop, the structures 1 and 2 should be
    !! identical to the input structures
    !!


    !!
    !! set found data
    !!
    !! center vectors
    call select_rc( nat1, coords1, c1min, rc1 )
    call select_rc( nat2, coords2, c2min, rc2 )
    !!
    !! rotation
    invb = transpose(beta_min)
    rotation = matmul(invb, gamma_min)
    !!
    !! translation
    translation(:) = rc1 - matmul(rotation, rc2)

    !! rotate + translate struc2
    do i = 1, nat2
      coords2(:,i) = matmul(rotation, coords2(:,i)) + translation
    end do

    !! permutation
    if( idxm .ne. 0 ) then
      !! if no bug in gamma (which happens if dist_k too small):
      !! find final permutations
      call cshda( nat1, typ1, coords1, &
            nat2, typ2, coords2, &
            999.9, found, dists )
      !!
      !! set final permutation
      !!
      permutation(:) = found(:)
    else
      !! all dists large
      dists(:) = 999.9
    endif

    !! set hd
    hd_out = maxval(dists(1:nat1))


  end subroutine ira_unify


  subroutine ira_svd( nat1, typ1_in, coords1_in, &
                      nat2, typ2_in, coords2_in, &
                      kmax_factor, rotation, translation, permutation, &
                      hd, rmsd )
    implicit none
    integer, intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1_in
    real, dimension(3,nat1), intent(in) :: coords1_in
    integer, intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2_in
    real, dimension(3,nat2), intent(in) :: coords2_in
    real, intent(in) :: kmax_factor
    real, dimension(3,3), intent(out) :: rotation
    real, dimension(3), intent(out) :: translation
    integer, dimension(nat2), intent(out) :: permutation
    real, intent(out) :: hd
    real, intent(out) :: rmsd

    integer, dimension(nat1) :: typ1
    real, dimension(3,nat1) :: coords1
    integer, dimension(nat2) :: typ2
    real, dimension(3,nat2) :: coords2
    integer, allocatable :: candidate_1(:), candidate_2(:)
    integer :: i
    real :: hd_out
    real, dimension(3,3) :: svd_rot
    real, dimension(3) :: rdum, svd_tr


    !! copy input data
    typ1(:) = typ1_in(:)
    typ2(:) = typ2_in(:)
    coords1(:,:) = coords1_in(:,:)
    coords2(:,:) = coords2_in(:,:)


    !!
    !! form candidates for central atms in 1 and 2
    !!
    allocate( candidate_1(1:nat1) )
    allocate( candidate_2(1:nat2) )
    call set_candidates( nat1, typ1, coords1, &
                        nat2, typ2, coords2, &
                        candidate_1, candidate_2 )

    !!
    !! call main ira with candidates
    !!
    call ira_unify( nat1, typ1, coords1, candidate_1, &
        nat2, typ2, coords2, candidate_2, &
        kmax_factor, rotation, translation, permutation, hd_out)

    !!
    !! apply found transformation
    !!
    typ2(:) = typ2(permutation(:))
    coords2(:,:) = coords2(:,permutation(:))
    !!
    do i = 1, nat2
      coords2(:,i) = matmul(rotation, coords2(:,i)) + translation
    end do

    !!
    !! call SVD
    !!
    call svdrot_m( nat1, typ1, coords1, &
        nat1, typ2(1:nat1), coords2(:,1:nat1), &
        svd_rot, svd_tr )

    !!
    !! apply svd
    !!
    ! do i = 1, nat2
    !    coords2(:,i) = matmul( svd_rot, coords2(:,i) ) + svd_tr
    ! end do

    !! rotate back to orig
    do i = 1, nat2
      coords2(:,i) = matmul(transpose(rotation),coords2(:,i)) - matmul(transpose(rotation),translation)
    end do

    !!
    !! put together apx and svd
    !!
    rotation = matmul(svd_rot, rotation)
    translation = matmul(svd_rot,translation) + svd_tr
    !! apply
    do i = 1, nat2
      coords2(:,i) = matmul(rotation,coords2(:,i)) + translation
    end do
    ! do i = 1, nat1
    !    coords1(:,i) = matmul(transpose(rotation),coords1(:,i)) - matmul(transpose(rotation),translation)
    ! end do

    !!
    !! compute dH and rmsd distances for output
    !!
    hd = 0.0
    do i = 1, nat1
      hd = max( hd, norm2(coords1(:,i) - coords2(:,i)) )
    end do

    rmsd = 0.0
    do i = 1, nat1
      rdum = coords1(:,i) - coords2(:,i)
      rmsd = rmsd + dot_product(rdum,rdum)
    end do
    rmsd = sqrt(rmsd/nat1)

  end subroutine ira_svd
