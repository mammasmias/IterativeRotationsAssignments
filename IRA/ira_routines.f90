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


  subroutine determinant( a, d )
    !! determinant of a 3x3 matrix a
    implicit none
    real, dimension(3,3), intent(in) :: a
    real, intent(out) :: d

    d = a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) + &
          a(1,2)*a(2,3)*a(3,1) - a(2,2)*a(3,1)*a(1,3) + &
          a(2,1)*a(3,2)*a(1,3) - a(3,3)*a(1,2)*a(2,1)

  end subroutine determinant


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


  subroutine cross_prod( a, b, c )
    !> @brief Cross product of two vectors
    implicit none
    real, dimension(3), intent(in) :: a
    real, dimension(3), intent(in) :: b
    real, dimension(3), intent(out) :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine cross_prod



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
    call cross_prod( basis(1,:), basis(2,:), basis(3,:) )
    norm_v = sqrt( dot_product(basis(3,:), basis(3,:)) )
    basis(3,:) = basis(3,:) / norm_v

    !!
    !! check for NaN. Can it happen?
    !!
    if( any(basis .ne. basis) ) then
       fail = .true.
       ! write(*,*) 'fail NaN'
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
                         kmax, gamma, gamma_idx, hd_out )
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
               nat2, typ2, coords2, found, dists )
          !!
          !! Hausdorff
          hd = maxval(dists(1:nat1) )
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
               nat2, typ2, coords2, found, dists )
          !!
          !! Hausdorff
          hd = maxval(dists(1:nat1))
          !!
          ! write(*,'(5i3,3f13.7)') i,j,nint(d_o(2,i)), nint(d_o(2,j)),m,hd, norm2(coords2(:,j))
          if( hd .lt. hd_old ) then
             hd_old = hd
             idx1 = nint(d_o(2,i))
             idx1 = i
             idx2 = nint(d_o(2,j))
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
    call set_orthonorm_bas( coords2(:,idx1), coords2(:,idx2), gamma, fail)
    gamma(3,:) = m_fin*gamma(3,:)

    !! for gamma_idx, use unsorted indices, which correspond to input order
    gamma_idx(1) = m_fin
    gamma_idx(2) = nint(d_o(2,idx1))
    gamma_idx(3) = nint(d_o(2,idx2))

    hd_out = hd_old
    deallocate( d_o )
    deallocate( found, dists )

    ! write(*,*) 'nbas tested:', count
    ! write(*,*) 'gamma idx:',gamma_idx(:)
    return
  end subroutine get_gamma_m


  subroutine ira_equal( nat1, typ1_in, coords1_in, &
       nat2, typ2_in, coords2_in, rotation, translation, permutation )

    !> @detail
    !! IRA + CShDA algorithm:
    !! Find approximate rotation and permutation by the basis rotation search,
    !! when the structures contain equal number of atoms.
    !!
    !! On output give rotation matrix, translation vector and permuttaion order,
    !! but don't actually rotate/translate/permute here!
    !!
    !! the output can be applied:
    !!
    !! matmul( rotation, coords2 ) + translation
    !! coords2(:,:) = coords2(:,permutation)
    !!
    implicit none
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1_in
    real, dimension(3,nat1),  intent(in) :: coords1_in
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2_in
    real, dimension(3,nat2),  intent(in) :: coords2_in
    real, dimension(3,3),     intent(out) :: rotation
    real, dimension(3),       intent(out) :: translation
    integer, dimension(nat1), intent(out) :: permutation
    !!
    integer :: i, j
    integer, allocatable :: typ1(:), typ2(:)
    real, allocatable :: coords1(:,:), coords2(:,:)
    real, dimension(3,3) :: gamma, beta, invb, rmat
    real, dimension(3) :: rc1, rc2
    real, allocatable :: dists(:)
    integer, allocatable :: found(:)
    real :: dist_k, hd_out
    logical :: fail
    real, allocatable :: d_o(:,:)
    integer, dimension(3) :: gamma_idx

    !! allocate working copies
    allocate( typ1(1:nat1), source = typ1_in)
    allocate( typ2(1:nat2), source = typ2_in)
    allocate( coords1(1:3,1:nat1), source = coords1_in)
    allocate( coords2(1:3,1:nat2), source = coords2_in)

    !!==========================================
    !! Set the centers as geometrical centers.
    !! If a better guess for a common point is present, modify this part.
    !! Keep in mind that rc1 and rc2 are needed for final translation vector.
    !!
    !! center of c1
    rc1=sum(coords1(:,:),2)/nat1
    do i = 1, nat1
       coords1(:,i) = coords1(:,i) - rc1
    end do

    !! center of c2
    rc2=sum(coords2(:,:),2)/nat2
    do i = 1, nat2
       coords2(:,i) = coords2(:,i) - rc2
    end do
    !!==========================================


    !! sort by size (distance from center)
    allocate( d_o(1:2,1:nat1))
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

    !! set initial permutation
    do i = 1, nat2
       permutation(i) = i
    end do

    deallocate( d_o )

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
          call set_orthonorm_bas( coords1(:,i),coords1(:,j), beta, fail)
          !! exit on first non-failed basis
          if( .not. fail ) exit
       end do
       if( .not. fail ) exit
    end do
    ! write(*,*) '1 in',i,j, norm2(coords1(:,i)), norm2(coords1(:,j))
    !!
    !! Measure "cutoff" of basis, multiply by a factor which is 1.2 by default.
    !! This factor should be greather than 1 in all cases. Decreasing it makes
    !! the search space for basis in structure 2 smaller, which can make the
    !! algorithm faster when structures are quite similar. It can also make it
    !! miss the optimal basis, in cases when the distortion is exactly on atoms
    !! corresponding to those which set the basis in structure 1.
    !!
    dist_k = max( norm2(coords1(:,i)), norm2(coords1(:,j)) )
    dist_k = dist_k*1.2
    !!
    !! rotate 1 temporarily to that bas
    !!
    do i = 1, nat1
       coords1(:,i) = matmul(beta, coords1(:,i))
    end do


    !!
    !! get basis in 2 (main IRA loop)
    !!
    allocate( found(1:nat1))
    allocate( dists(1:nat1))
    call get_gamma_m(nat1, typ1, coords1, &
                     nat2, typ2, coords2, &
                     dist_k, gamma, gamma_idx, hd_out )
    !!
    !! rotate to found basis transpose(gamma), coords1 is already in beta
    !!
    do i = 1, nat1
       coords1(:,i) = matmul(transpose(gamma),coords1(:,i))
    end do
    !!
    !! find permutations (this could be skipped)
    !!
    call cshda( nat1, typ1, coords1, &
         nat2, typ2, coords2, found, dists )
    !!
    !! permute
    !!
    ! typ2(:) = typ2(found(1:nat2))
    ! coords2(:,:) = coords2(:,found(1:nat2))
    call permute_int_1d(nat2, typ2, found(1:nat2) )
    call permute_real_2d(nat2, 3, coords2, found(1:nat2) )

    ! write(*,*) nat1
    ! write(*,*) 'c1'
    ! do k = 1, nat1
    !    write(*,*) typ1(k), coords1(:,k)
    ! end do
    ! write(*,*) nat2
    ! write(*,*) 'c2'
    ! do k = 1, nat2
    !    write(*,*) typ2(k), coords2(:,k)
    ! end do
    !!
    !! set final R_apx matrix
    !!
    invb=transpose(beta)
    rmat = matmul(invb,gamma)
    !!
    !! approx rotation and translation
    !!
    rotation(:,:) = rmat(:,:)
    translation(:) = rc1-matmul(rmat,rc2)

    !!
    !! re-set originals, to check if all is ok
    !!
    typ1 = typ1_in
    typ2 = typ2_in
    coords1 = coords1_in
    coords2 = coords2_in

    !! rotate+translate
    do i = 1, nat2
       coords2(:,i) = matmul(rotation,coords2(:,i)) + translation
    end do

    !! find final permutations
    call cshda( nat1, typ1, coords1, &
         nat2, typ2, coords2, found, dists )
    !!
    !! set final permutation
    !!
    permutation(:) = found(:)


    deallocate( coords1, coords2 )
  end subroutine ira_equal



  subroutine ira_nonequal( nat1, typ1_in, coords1_in, &
                           nat2, typ2_in, coords2_in, &
                           c_atm1, rotation, translation, permutation )

    !> @detail
    !! IRA + CShDA algorithm:
    !! Find approximate rotation and permutation by the basis rotation search,
    !! for each candidate central atom in structure 2.
    !!
    !! On output give rotation matrix, translation vector and permuttaion order,
    !! but don't actually rotate/translate/permute here!
    !!
    !! the output can be applied:
    !!
    !! matmul( rotation, coords2 ) + translation
    !! coords2(:,:) = coords2(:,permutation)
    !!
    implicit none
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1_in
    real, dimension(3,nat1),  intent(in) :: coords1_in
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2_in
    real, dimension(3,nat2),  intent(in) :: coords2_in
    integer,                  intent(in) :: c_atm1
    real, dimension(3,3),     intent(out) :: rotation
    real, dimension(3),       intent(out) :: translation
    integer, dimension(nat2), intent(out) :: permutation
    !!
    integer :: i, j
    integer, allocatable :: typ1(:), typ2(:)
    real, allocatable :: coords1(:,:), coords2(:,:)
    real, dimension(3,3) :: gamma, beta, invb, rmat
    real, dimension(3) :: rc1, rc2
    real, allocatable :: dists(:)
    integer, allocatable :: found(:)
    real :: hd, hd_old, dist_k
    integer :: idx1, idx2, ij
    logical :: fail
    real, allocatable :: d_o(:,:)
    integer :: tc, idxm, typ_c1
    integer, dimension(3) :: gamma_idx

    !! allocate working copies
    allocate( typ1(1:nat1), source = typ1_in)
    allocate( typ2(1:nat2), source = typ2_in)
    allocate( coords1(1:3,1:nat1), source = coords1_in)
    allocate( coords2(1:3,1:nat2), source = coords2_in)

    !! center of c1 is atom index c_atm1 from input
    rc1=coords1(:,c_atm1)
    typ_c1 = typ1(c_atm1)
    do i = 1, nat1
       coords1(:,i) = coords1(:,i) - rc1
    end do

    !! sort by size
    allocate( d_o(1:2,1:nat1))
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

    !! set initial permutation
    do i = 1, nat2
       permutation(i) = i
    end do

    deallocate( d_o )

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
          call set_orthonorm_bas( coords1(:,i),coords1(:,j), beta, fail)
          !! exit on first non-failed basis
          if( .not. fail ) exit
       end do
       if( .not. fail ) exit
    end do
    ! write(*,*) '1 in',i,j, norm2(coords1(:,i)), norm2(coords1(:,j))
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
    dist_k = dist_k*1.2
    !!
    !! rotate 1 to that bas
    !!
    do i = 1, nat1
       coords1(:,i) = matmul(beta, coords1(:,i))
    end do
    !!
    !! get basis in 2 (main IRA loop), over trial central atoms
    !!
    allocate( found(1:nat2))
    allocate( dists(1:nat2))
    hd_old = 9999.9
    !!
    do ij = 1, nat2
       !!
       !! both central atoms should be of same typ
       if( typ2(ij) .ne. typ_c1 ) cycle
       !!
       !! shift to trial central atom
       !!
       rc2 = coords2(:,ij)
       !!
       do i = 1, nat2
          coords2(:,i) = coords2(:,i) - rc2
       end do
       !!
       !! get gamma for this central atm
       !!
       call get_gamma_m(nat1, typ1(1:nat1), coords1(1:3,1:nat1), &
                        nat2, typ2(1:nat2), coords2(1:3,1:nat2), &
                        dist_k, gamma, gamma_idx, hd )
       !!
       !! rotate to found basis
       !!
       do i = 1, nat1
          coords1(:,i) = matmul(transpose(gamma),coords1(:,i))
       end do

       !!
       !! find final permutations
       !!
       call cshda( nat1, typ1(1:nat1), coords1(1:3,1:nat1), &
                   nat2, typ2(1:nat2), coords2(1:3,1:nat2), &
                   found(1:nat2), dists(1:nat2) )
       !!
       !! Hausdorff of this central, up to nat1
       hd = maxval(dists(1:nat1))
       !!
       !! store minimum
       if( hd .lt. hd_old ) then
          hd_old = hd
          tc = ij
          idxm = gamma_idx(1)
          idx1 = gamma_idx(2)
          idx2 = gamma_idx(3)
       endif

       ! write(*,*) ij, hd

       !! rotate back to orig
       do i = 1, nat1
          coords1(:,i) = matmul(gamma, coords1(:,i))
       end do

       !! shift back to orig
       do i = 1, nat2
          coords2(:,i) = coords2(:,i) + rc2
       end do

    end do
    ! write(*,*) hd_old, tc, idxm, idx1, idx2


    !! translate to found central atom
    rc2 = coords2(:,tc)
    do i = 1, nat2
       coords2(:,i) = coords2(:,i) - rc2
    end do

    !! set the found gamma matrix
    call set_orthonorm_bas( coords2(:,idx1),coords2(:,idx2), gamma, fail)
    gamma(3,:) = gamma(3,:)*idxm

    !! rotate to transpose(gamma), coords1 is already in beta.
    do i = 1, nat1
       coords1(:,i) = matmul(transpose(gamma),coords1(:,i))
    end do

    !! get permutation (this could be skipped)
    call cshda( nat1, typ1(1:nat1), coords1(1:3,1:nat1), &
         nat2, typ2(1:nat2), coords2(1:3,1:nat2), found(1:nat2), dists(1:nat2) )
    !! permute
    ! typ2(:) = typ2(found(1:nat2))
    ! coords2(:,:) = coords2(:,found(1:nat2))
    call permute_int_1d(nat2, typ2, found(1:nat2) )
    call permute_real_2d(nat2, 3, coords2, found(1:nat2) )


    ! write(*,*) nat1
    ! write(*,*) 'c1'
    ! do k = 1, nat1
    !    write(*,*) typ1(k), coords1(:,k)
    ! end do
    ! write(*,*) nat2
    ! write(*,*) 'c2'
    ! do k = 1, nat2
    !    write(*,*) typ2(k), coords2(:,k)
    ! end do
    !!
    !! set final R_apx matrix
    !!
    invb=transpose(beta)
    rmat = matmul(invb,gamma)
    !!
    !! approx rotation and translation
    !!
    rotation(:,:) = rmat(:,:)
    translation(:) = rc1-matmul(rmat,rc2)

    !!
    !! re-set originals, to check if all is ok
    !!
    typ1 = typ1_in
    typ2 = typ2_in
    coords1 = coords1_in
    coords2 = coords2_in

    !! rotate+translate
    do i = 1, nat2
       coords2(:,i) = matmul(rotation,coords2(:,i)) + translation
    end do

    !! find final permutations
    call cshda( nat1, typ1, coords1, &
         nat2, typ2, coords2, found, dists )
    ! write(*,*) 'found'
    ! do i = 1, nat2
    !    write(*,*) i, found(i), dists(i)
    ! end do

    !!
    !! set final permutation
    !!
    permutation(:) = found(:)

    deallocate( coords1, coords2 )
  end subroutine ira_nonequal

