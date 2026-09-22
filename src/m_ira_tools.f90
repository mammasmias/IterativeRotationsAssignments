
!> @cond SKIP
module m_ira_tools
  use ira_precision
  implicit none
  public

contains

  !> @brief determinant of a 3x3 matrix a
  pure subroutine determinant3x3(a, d)
    use ira_precision
    implicit none
    real(rp), dimension(3, 3), intent(in) :: a
    real(rp), intent(out) :: d

    d = a(1, 1)*a(2, 2)*a(3, 3) - a(1, 1)*a(2, 3)*a(3, 2) + &
         a(1, 2)*a(2, 3)*a(3, 1) - a(2, 2)*a(3, 1)*a(1, 3) + &
         a(2, 1)*a(3, 2)*a(1, 3) - a(3, 3)*a(1, 2)*a(2, 1)

  end subroutine determinant3x3

  subroutine periodic(c)
    !--------------------------------
    ! periodic boundary condition, for 3 dimensional vector input in crist coords.
    !--------------------------------
    implicit none
    real(rp), dimension(3),intent(inout) :: c
    integer(ip) :: i

    do i = 1, 3
       if( c(i) < -0.5_rp ) c(i) = c(i) + 1.0_rp
       if( c(i) >= 0.5_rp ) c(i) = c(i) - 1.0_rp
    end do

  end subroutine periodic


  subroutine cart_to_crist(xpp,ct)
    !!----------------------------
    !! cartesian to crystallographic coordinates transform, in 3-dimension
    !! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
    !! This routine does the transpose of B implicitly
    !! --------
    !! xpp(3)      ==> input vector of position in cartesian
    !! ct(3,3)     ==> conversion matrix, vectors of the Bravais lattice in rows
    !!
    !!      ct = a1 a2 a3
    !!           b1 b2 b3
    !!           c1 c2 c3
    !!----------------------------
    !! bt(3,3) ==> inverse matrix of ct, used locally
    !! xc(3)   ==> copy of xpp, used locally
    !! detct   ==> determinant of ct, used locally
    !!
    implicit none
    real(rp), dimension(3),   intent(inout) :: xpp
    real(rp), dimension(3,3), intent(in)    :: ct

    real(rp),dimension(3) :: xc
    real(rp) :: detct
    real(rp), dimension(3,3) :: bt

    ! -----------------------------------------------
    !  inverse matrix of ct(:,:)
    !------------------------------------------------
    detct=ct(1,1)*ct(2,2)*ct(3,3)+&
         ct(1,2)*ct(2,3)*ct(3,1)+&
         ct(2,1)*ct(3,2)*ct(1,3)&
         -ct(1,3)*ct(2,2)*ct(3,1)&
         -ct(3,2)*ct(2,3)*ct(1,1)&
         -ct(1,2)*ct(2,1)*ct(3,3)

    bt(1,1)= ct(2,2)*ct(3,3)-ct(2,3)*ct(3,2)
    bt(1,2)=-(ct(1,2)*ct(3,3)-ct(1,3)*ct(3,2))
    bt(1,3)= ct(1,2)*ct(2,3)-ct(1,3)*ct(2,2)
    bt(2,1)=-(ct(2,1)*ct(3,3)-ct(2,3)*ct(3,1))
    bt(2,2)= ct(1,1)*ct(3,3)-ct(3,1)*ct(1,3)
    bt(2,3)=-(ct(1,1)*ct(2,3)-ct(1,3)*ct(2,1))
    bt(3,1)= ct(2,1)*ct(3,2)-ct(2,2)*ct(3,1)
    bt(3,2)=-(ct(1,1)*ct(3,2)-ct(1,2)*ct(3,1))
    bt(3,3)= ct(1,1)*ct(2,2)-ct(2,1)*ct(1,2)
    !------------------------------------------------

    xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))/detct
    xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))/detct
    xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))/detct

    xpp(:) = xc(:)

  end subroutine cart_to_crist


  subroutine crist_to_cart(xpp,bt)
    !!--------------------------------
    !! crystallographic to cartesian transformation in 3-dimensions
    !! R_cart = B * v_crist; where B is the matrix formed by cell vectors vertically
    !! This routine does the transpose implicitly!
    !! -----------
    !! xpp(3)    ==> input vector in crystallographic, output vector in cartesian
    !! bt(3,3)   ==> input conversion matrix, vectors of the Bravais lattice in rows
    !!
    !!         bt = a1 a2 a3
    !!              b1 b2 b3
    !!              c1 c2 c3
    !!-----
    !! xc(3)   ==> local vector
    !!
    implicit none

    real(rp), dimension(3),   intent(inout) :: xpp

    real(rp), dimension(3,3), intent(in)    :: bt
    real(rp), dimension(3) :: xc


    xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))
    xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))
    xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))

    xpp(:) = xc(:)

  end subroutine crist_to_cart


  subroutine pbc_vec( vec, lat )
    !! apply pbc of lattice 'lat' to a vector 'vec'
    implicit none
    real(rp), dimension(3), intent(inout) :: vec
    real(rp), dimension(3,3), intent(in) :: lat

    call cart_to_crist( vec, lat )
    call periodic( vec )
    call crist_to_cart( vec, lat )

  end subroutine pbc_vec


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
    use m_ira_error
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



end module m_ira_tools
!> @endcond
