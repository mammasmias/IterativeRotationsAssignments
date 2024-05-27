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


module sofi_tools

  implicit none
  public

  !! maximum size of output lists. The real number of elements is output as integer
  !! This is pre-defined for security.
  integer, parameter :: nmax = 200


  !! matrix-distance thr: value 1.4 captures 1/6*2pi rotations, which are
  !! needed to distinguish D2h from D6h, since D2h is subgroup of D6h!
  !! see the comment under matrix_distance in sofi_tools.f90 for some ref. values
  !! 1.07 captures 1/8*2pi
  !! NOTE: decreasing this value gives "higher resolution", but slows down the algo!
  !!   Also, groups with order > 6 are super rare in atomic clusters. But can happen in
  !!   for example nanotubes, where main ax is in center of tube, around this ax
  !!   many rotations can happen, then order of group can be any.
  ! m_thr = 1.4
  ! real, parameter :: m_thr = 1.07
  real, parameter :: m_thr = 0.73
  ! real, parameter :: m_thr = 0.58
  ! real, parameter :: m_thr = 0.35


  !! Schoenflies symbols for operations
  character(len=1), parameter :: &
       OP_ERROR      = "X", &
       OP_IDENTITY   = "E", &
       OP_INVERSION  = "I", &
       OP_PROP_ROT   = "C", &
       OP_IMPROP_ROT = "S", &
       OP_MIRROR     = "S" !! "M"


  real, parameter :: pi = 4.0*atan(1.0)
  real, parameter :: epsilon = 1e-6
  real, parameter :: collinearity_thr = 0.95

contains

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


  recursive function gcd_rec(u, v) result(gcd)
    !! greatest common denominator
    integer             :: gcd
    integer, intent(in) :: u, v

    if (mod(u, v) /= 0) then
       gcd = gcd_rec(v, mod(u, v))
    else
       gcd = v
    end if
  end function gcd_rec


  subroutine matrix_distance( a, b, dist )
    !! distance between two 3x3 matrices a and b, computed element-wise as sqrt( sum( (a-b)^2 ) )
    !! This is an arbitrary choice.
    !!
    !! Some reference values of distance:
    !!
    !!  a, b rotated by angle    ||   dist
    !!               1/2 * 2pi   ||   2.8
    !!               1/3 * 2pi   ||   2.45
    !!               1/6 * 2pi   ||   1.4142
    !!              1/12 * 2pi   ||   0.7321
    !!              1/18 * 2pi   ||   0.4912
    !!              1/24 * 2pi   ||   0.3692
    implicit none
    real, dimension(3,3), intent(in) :: a, b
    real, intent(out) :: dist

    integer :: i, j

    dist = 0.0
    do i = 1, 3
       do j  =1, 3
          dist = dist + ( a(i,j) - b(i,j) )**2
       end do
    end do
    dist = sqrt(dist)
    ! call md(a,b,dist)

  end subroutine matrix_distance

  subroutine md(a,b,dist)
    !! unused
    implicit none
    real, dimension(3,3), intent(in) :: a, b
    real, intent(out) :: dist

    real, dimension(3,3) :: rij
    real :: tr

    rij = matmul( a, transpose(b) )
    rij = transpose(rij)
    tr = rij(1,1) + rij(2,2) + rij(3,3) - 3.0
    tr = abs(tr)
    dist = tr

  end subroutine md


  SUBROUTINE diag(n, A, eigvals, vec)
    !> @brief
    !! assuming a general square matrix (can be nonsymmetric).
    !! On output A is overwritten by eigenvectors in rows, if vec=0, then
    !! A is just 0.0 on output.
    !!
    !! @param [in]    n	      dimension of matrix A
    !! @param [inout] A	      matrix to be diagonalised, overwritten by eigenvectors in columns on output
    !! @param [out]   eigvals   output vector of eigenvalues, not sorted!
    !! @param [in]    vec	      0 if don't want to compute eigenvectors, 1 otherwise
    !!
    IMPLICIT NONE
    INTEGER,              intent(in) :: n
    REAL, DIMENSION(n,n), intent(inout) :: A
    REAL, DIMENSION(n),   intent(out) :: eigvals
    INTEGER,              intent(in) :: vec
    REAL, DIMENSION(n) :: eigvals_i !! imaginary part of the eigenvalues
    REAL, DIMENSION(n,n) :: eigvec
    INTEGER :: lda
    INTEGER :: lwork
    REAL :: Dummy(1000)
    INTEGER :: info
    CHARACTER(len=1) :: getvec

    real, dimension(3,3) :: mat
    getvec = 'N'
    if( vec == 1 ) getvec='V'
    lda = n
    eigvals_i(:) = 0.0
    eigvec(:,:) = 0.0

    mat = A
    !! test workspace
    lwork = -1
    call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
              dummy, 1, eigvec, n, dummy, lwork, info)

    !! choose optimal size of workspace (as in example from intel website)
    lwork = min( 1000, nint(dummy(1)) )
    !! compute stuffs
    call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
                dummy, 1, eigvec, n, dummy, lwork, info)

    ! write(*,*) 'diagonal matrix'
    ! ! mat = matmul(mat, transpose(A))
    ! write(*,'(3f9.4)')a(1,:)
    ! write(*,'(3f9.4)')a(2,:)
    ! write(*,'(3f9.4)')a(3,:)

    ! write(*,*) 'eigvals'
    ! write(*,'(3f9.5)') eigvals
    ! write(*,'(3f9.5)') eigvals_i
    ! write(*,*) eigvec(:,1)
    ! write(*,*) matmul(mat,eigvec(:,1))
    ! write(*,*)
    ! write(*,*) eigvec(:,2)
    ! write(*,*) matmul(mat,eigvec(:,2))
    ! write(*,*)
    ! write(*,*) eigvec(:,3)
    ! write(*,*) matmul(mat,eigvec(:,3))

    !! overwrite a on output
    A(:,:) = eigvec(:,:)
  END SUBROUTINE diag


  function check_pg_Nop( pg, Nop ) result( err )
    !! check whether the name of PG is consistent with the number of operations Nop of that group
    !! NOTE: errors currently go up to order 10 PGs.
    implicit none
    integer :: err
    character(len=10), intent(in) :: pg
    integer, intent(in) :: Nop

    err = 0
    !! if Nop is less than should be, return err=-1
    if( &
         pg == 'C1'   .and. Nop .lt. 1   .or. &
         pg == 'Cs'   .and. Nop .lt. 2   .or. &
         pg == 'C2'   .and. Nop .lt. 2   .or. &
         pg == 'C2h'  .and. Nop .lt. 4   .or. &
         pg == 'C2v'  .and. Nop .lt. 4   .or. &
         pg == 'C3'   .and. Nop .lt. 3   .or. &
         pg == 'C3h'  .and. Nop .lt. 6   .or. &
         pg == 'C3v'  .and. Nop .lt. 6   .or. &
         pg == 'C4'   .and. Nop .lt. 4   .or. &
         pg == 'C4h'  .and. Nop .lt. 8   .or. &
         pg == 'C4v'  .and. Nop .lt. 8   .or. &
         pg == 'C5'   .and. Nop .lt. 5   .or. &
         pg == 'C5h'  .and. Nop .lt. 10  .or. &
         pg == 'C5v'  .and. Nop .lt. 10  .or. &
         pg == 'C6'   .and. Nop .lt. 6   .or. &
         pg == 'C6h'  .and. Nop .lt. 12  .or. &
         pg == 'C6v'  .and. Nop .lt. 12  .or. &
         pg == 'C7'   .and. Nop .lt. 7   .or. &
         pg == 'C7h'  .and. Nop .lt. 14  .or. &
         pg == 'C7v'  .and. Nop .lt. 14  .or. &
         pg == 'C8'   .and. Nop .lt. 8   .or. &
         pg == 'C8h'  .and. Nop .lt. 16  .or. &
         pg == 'C8v'  .and. Nop .lt. 16  .or. &
         pg == 'C9'   .and. Nop .lt. 9   .or. &
         pg == 'C9h'  .and. Nop .lt. 18  .or. &
         pg == 'C9v'  .and. Nop .lt. 18  .or. &
         pg == 'C10'  .and. Nop .lt. 10  .or. &
         pg == 'C10h' .and. Nop .lt. 20  .or. &
         pg == 'C10v' .and. Nop .lt. 20  .or. &
         pg == 'D2'   .and. Nop .lt. 4   .or. &
         pg == 'D2h'  .and. Nop .lt. 8   .or. &
         pg == 'D2d'  .and. Nop .lt. 8   .or. &
         pg == 'D3'   .and. Nop .lt. 6   .or. &
         pg == 'D3h'  .and. Nop .lt. 12  .or. &
         pg == 'D3d'  .and. Nop .lt. 12  .or. &
         pg == 'D4'   .and. Nop .lt. 8   .or. &
         pg == 'D4h'  .and. Nop .lt. 16  .or. &
         pg == 'D4d'  .and. Nop .lt. 16  .or. &
         pg == 'D5'   .and. Nop .lt. 10  .or. &
         pg == 'D5h'  .and. Nop .lt. 20  .or. &
         pg == 'D5d'  .and. Nop .lt. 20  .or. &
         pg == 'D6'   .and. Nop .lt. 12  .or. &
         pg == 'D6h'  .and. Nop .lt. 24  .or. &
         pg == 'D6d'  .and. Nop .lt. 24  .or. &
         pg == 'D7'   .and. Nop .lt. 14  .or. &
         pg == 'D7h'  .and. Nop .lt. 28  .or. &
         pg == 'D7d'  .and. Nop .lt. 28  .or. &
         pg == 'D8'   .and. Nop .lt. 16  .or. &
         pg == 'D8h'  .and. Nop .lt. 32  .or. &
         pg == 'D8d'  .and. Nop .lt. 32  .or. &
         pg == 'D9'   .and. Nop .lt. 18  .or. &
         pg == 'D9h'  .and. Nop .lt. 36  .or. &
         pg == 'D9d'  .and. Nop .lt. 36  .or. &
         pg == 'D10'  .and. Nop .lt. 20  .or. &
         pg == 'D10h' .and. Nop .lt. 40  .or. &
         pg == 'D10d' .and. Nop .lt. 40  .or. &
         pg == 'T'    .and. Nop .lt. 12  .or. &
         pg == 'Th'   .and. Nop .lt. 24  .or. &
         pg == 'Td'   .and. Nop .lt. 24  .or. &
         pg == 'O'    .and. Nop .lt. 24  .or. &
         pg == 'Oh'   .and. Nop .lt. 48  .or. &
         pg == 'I'    .and. Nop .lt. 60  .or. &
         pg == 'Ih'   .and. Nop .lt. 120 .or. &
         pg == 'Ci'   .and. Nop .lt. 2   .or. &
         pg == 'S4'   .and. Nop .lt. 4   .or. &
         pg == 'S6'   .and. Nop .lt. 6   .or. &
         pg == 'S8'   .and. Nop .lt. 8   .or. &
         pg == 'S10'  .and. Nop .lt. 10   ) err = -1

    !! if Nop is more than should be, return err=1
    if( &
         pg == 'C1'   .and. Nop .gt. 1   .or. &
         pg == 'Cs'   .and. Nop .gt. 2   .or. &
         pg == 'C2'   .and. Nop .gt. 2   .or. &
         pg == 'C2h'  .and. Nop .gt. 4   .or. &
         pg == 'C2v'  .and. Nop .gt. 4   .or. &
         pg == 'C3'   .and. Nop .gt. 3   .or. &
         pg == 'C3h'  .and. Nop .gt. 6   .or. &
         pg == 'C3v'  .and. Nop .gt. 6   .or. &
         pg == 'C4'   .and. Nop .gt. 4   .or. &
         pg == 'C4h'  .and. Nop .gt. 8   .or. &
         pg == 'C4v'  .and. Nop .gt. 8   .or. &
         pg == 'C5'   .and. Nop .gt. 5   .or. &
         pg == 'C5h'  .and. Nop .gt. 10  .or. &
         pg == 'C5v'  .and. Nop .gt. 10  .or. &
         pg == 'C6'   .and. Nop .gt. 6   .or. &
         pg == 'C6h'  .and. Nop .gt. 12  .or. &
         pg == 'C6v'  .and. Nop .gt. 12  .or. &
         pg == 'C7'   .and. Nop .gt. 7   .or. &
         pg == 'C7h'  .and. Nop .gt. 14  .or. &
         pg == 'C7v'  .and. Nop .gt. 14  .or. &
         pg == 'C8'   .and. Nop .gt. 8   .or. &
         pg == 'C8h'  .and. Nop .gt. 16  .or. &
         pg == 'C8v'  .and. Nop .gt. 16  .or. &
         pg == 'C9'   .and. Nop .gt. 9   .or. &
         pg == 'C9h'  .and. Nop .gt. 18  .or. &
         pg == 'C9v'  .and. Nop .gt. 18  .or. &
         pg == 'C10'  .and. Nop .gt. 10  .or. &
         pg == 'C10h' .and. Nop .gt. 20  .or. &
         pg == 'C10v' .and. Nop .gt. 20  .or. &
         pg == 'D2'   .and. Nop .gt. 4   .or. &
         pg == 'D2h'  .and. Nop .gt. 8   .or. &
         pg == 'D2d'  .and. Nop .gt. 8   .or. &
         pg == 'D3'   .and. Nop .gt. 6   .or. &
         pg == 'D3h'  .and. Nop .gt. 12  .or. &
         pg == 'D3d'  .and. Nop .gt. 12  .or. &
         pg == 'D4'   .and. Nop .gt. 8   .or. &
         pg == 'D4h'  .and. Nop .gt. 16  .or. &
         pg == 'D4d'  .and. Nop .gt. 16  .or. &
         pg == 'D5'   .and. Nop .gt. 10  .or. &
         pg == 'D5h'  .and. Nop .gt. 20  .or. &
         pg == 'D5d'  .and. Nop .gt. 20  .or. &
         pg == 'D6'   .and. Nop .gt. 12  .or. &
         pg == 'D6h'  .and. Nop .gt. 24  .or. &
         pg == 'D6d'  .and. Nop .gt. 24  .or. &
         pg == 'D7'   .and. Nop .gt. 14  .or. &
         pg == 'D7h'  .and. Nop .gt. 28  .or. &
         pg == 'D7d'  .and. Nop .gt. 28  .or. &
         pg == 'D8'   .and. Nop .gt. 16  .or. &
         pg == 'D8h'  .and. Nop .gt. 32  .or. &
         pg == 'D8d'  .and. Nop .gt. 32  .or. &
         pg == 'D9'   .and. Nop .gt. 18  .or. &
         pg == 'D9h'  .and. Nop .gt. 36  .or. &
         pg == 'D9d'  .and. Nop .gt. 36  .or. &
         pg == 'D10'  .and. Nop .gt. 20  .or. &
         pg == 'D10h' .and. Nop .gt. 40  .or. &
         pg == 'D10d' .and. Nop .gt. 40  .or. &
         pg == 'T'    .and. Nop .gt. 12  .or. &
         pg == 'Th'   .and. Nop .gt. 24  .or. &
         pg == 'Td'   .and. Nop .gt. 24  .or. &
         pg == 'O'    .and. Nop .gt. 24  .or. &
         pg == 'Oh'   .and. Nop .gt. 48  .or. &
         pg == 'I'    .and. Nop .gt. 60  .or. &
         pg == 'Ih'   .and. Nop .gt. 120 .or. &
         pg == 'Ci'   .and. Nop .gt. 2   .or. &
         pg == 'S4'   .and. Nop .gt. 4   .or. &
         pg == 'S6'   .and. Nop .gt. 6   .or. &
         pg == 'S8'   .and. Nop .gt. 8   .or. &
         pg == 'S10'  .and. Nop .gt. 10   ) err = 1

    return

  end function check_pg_Nop


  function find_inversion( nbas, op ) result( has_inversion )
    !! find if there is inversion in the list
    implicit none
    logical :: has_inversion
    integer, intent(in) :: nbas
    character(len=1), dimension(nbas), intent(in) :: op
    integer :: i

    has_inversion = .false.
    do i = 1, nbas
       if( op(i) == OP_INVERSION ) then
          has_inversion = .true.
       end if
    end do
    return

  end function find_inversion

  function find_sigma( nbas, op, n_int ) result( has_sigma )
    !! find if there is sigma op in the list
    implicit none
    logical :: has_sigma
    integer, intent(in) :: nbas
    character(len=1), dimension(nbas), intent(in) :: op
    integer, dimension(nbas), intent(in) :: n_int
    integer :: i

    has_sigma = .false.
    do i = 1, nbas
       if( op(i) .eq. OP_IMPROP_ROT .and. n_int(i) .eq. 0 ) has_sigma = .true.
    end do
    return
  end function find_sigma

  function find_cn( nbas, op, n_int ) result( has_cn )
    !! find if there is Cn op in the list
    implicit none
    logical :: has_cn
    integer, intent(in) :: nbas
    character(len=1), dimension(nbas), intent(in) :: op
    integer, dimension(nbas), intent(in) :: n_int
    integer :: i

    has_cn = .false.
    do i = 1, nbas
       if( op(i) .eq. OP_PROP_ROT .and. n_int(i) .gt. 1 ) has_cn = .true.
    end do
    return

  end function find_cn


  function op_valid_ext_b(rmat, bfield) result( is_valid )
    !! test if Op matrix is valid in externam Bfield:
    !!
    !!  det(M) MB = B
    !!
    implicit none
    logical :: is_valid
    real, dimension(3,3), intent(in) :: rmat
    real, dimension(3), intent(in) :: bfield

    real, dimension(3) :: bdir, vec
    real :: det, dotp

    is_valid = .false.

    !! direction of B field
    bdir = bfield/norm2(bfield)

    call determinant3x3(rmat, det)

    !! det(M)*MB
    vec = matmul( rmat, bdir )
    vec = vec * det

    !! compare vec to bdir, shoudl be equal: dot = 1.0
    dotp = dot_product( vec, bdir )
    if( dotp .gt. 0.999 ) is_valid = .true.

  end function op_valid_ext_b

  subroutine construct_reflection( ax, rmat )
    !! construct reflection matrix from axis, as:
    !! rmat = I - 2.0* ( <x x^T> / <x^T x>), where x is the axis
    real, dimension(3), intent(in) :: ax
    real, dimension(3,3), intent(out) :: rmat

    real, dimension(3,3) :: id, rr
    integer :: i, j

    id(:,:) = 0.0
    do i = 1, 3
       id(i,i) = 1.0
    end do

    !! outer product
    forall (i=1:3)
       forall(j=1:3) rr(i,j) = ax(i)*ax(j)
    end forall
    rr(:,:) = rr / norm2(ax)

    rmat(:,:) = id(:,:) - 2.0*rr(:,:)

  end subroutine construct_reflection
  subroutine construct_rotation(ax_in, angle, rmat)
    !! wikipedia/rotation_matrix#quaternion
    implicit none
    real, dimension(3), intent(in) :: ax_in
    real, intent(in) :: angle
    real, dimension(3,3), intent(out) :: rmat

    real, dimension(3) :: ax
    real :: c, s, c_one

    ax = ax_in/norm2(ax_in)

    c = cos(angle)
    s = sin(angle)
    c_one = 1.0 - c

    rmat(1,1) = ax(1)**2*c_one + c
    rmat(1,2) = ax(1)*ax(2)*c_one - ax(3)*s
    rmat(1,3) = ax(1)*ax(3)*c_one + ax(2)*s

    rmat(2,1) = ax(2)*ax(1)*c_one + ax(3)*s
    rmat(2,2) = ax(2)**2*c_one + c
    rmat(2,3) = ax(2)*ax(3)*c_one - ax(1)*s

    rmat(3,1) = ax(3)*ax(1)*c_one - ax(2)*s
    rmat(3,2) = ax(3)*ax(2)*c_one + ax(1)*s
    rmat(3,3) = ax(3)**2*c_one + c

    !call determinant(rmat,dr)
    !write(*,*) 'det r:',dr
  end subroutine construct_rotation


  !> @details
  !! convention for axis direction:
  !! flip such that z>0
  !! if z==0, then flip such that x>0
  !! if x==0, then flip such that y>0
  subroutine ax_convention( ax )
    implicit none
    real, intent(inout) :: ax(3)

    real :: flip

    flip = 1.0
    if( ax(3) .lt. -epsilon ) then
       !! z is negative, flip
       flip = -flip
    elseif( abs(ax(3)) < epsilon ) then
       !! z==0, check x
       if( ax(1) < -epsilon ) then
          !! x is negative, flip
          flip = -flip
       elseif( abs(ax(1)) < epsilon ) then
          !! x==0, check y
          if( ax(2) < -epsilon ) then
             !! y is negative, flip
             flip = -flip
          end if
       end if
    end if
    ax = ax*flip
  end subroutine ax_convention


end module sofi_tools
