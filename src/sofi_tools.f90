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

  use ira_precision
  implicit none
  public

  !! maximum size of output lists. The real number of elements is output as integer
  !! This is pre-defined for security.
  integer(ip), parameter :: nmax = 400


  !! matrix-distance thr: value 1.4 captures 1/6*2pi rotations, which are
  !! needed to distinguish D2h from D6h, since D2h is subgroup of D6h!
  !! see the comment under matrix_distance in sofi_tools.f90 for some ref. values
  !! 1.07 captures 1/8*2pi
  !! NOTE: decreasing this value gives "higher resolution", but slows down the algo!
  !!   Also, groups with order > 8 are super rare in atomic clusters. But can happen in
  !!   for example nanotubes, where main ax is in center of tube, around this ax
  !!   many rotations can happen, then order of group can be any.
  ! real(rp), parameter :: m_thr = 1.4_rp      !! C6
  ! real(rp), parameter :: m_thr = 1.07_rp     !! C8
  ! real(rp), parameter :: m_thr = 0.73_rp     !! C12
  ! real(rp), parameter :: m_thr = 0.49_rp     !! C18
  ! real(rp), parameter :: m_thr = 0.36_rp     !! C24
  ! real(rp), parameter :: m_thr = 0.19_rp     !! C48
  ! real(rp), parameter :: m_thr = 0.092_rp     !! C96
  real(rp), parameter :: m_thr = 0.044_rp     !! C200


  !! Schoenflies symbols for operations
  character(len=1), parameter :: &
       OP_ERROR      = "X", &
       OP_IDENTITY   = "E", &
       OP_INVERSION  = "I", &
       OP_PROP_ROT   = "C", &
       OP_IMPROP_ROT = "S", &
       OP_MIRROR     = "S" !! "M"


  real(rp), parameter :: pi = 4.0_rp*atan(1.0_rp)
  real(rp), parameter :: epsilon = 1e-6_rp
  real(rp), parameter :: collinearity_thr = 0.95_rp

  !! limit value for n in sofi_analmat.
  ! integer(ip), parameter :: lim_n_val = 24
  ! integer(ip), parameter :: lim_n_val = 48
  ! integer(ip), parameter :: lim_n_val = 96
  integer(ip), parameter :: lim_n_val = 200


contains

  subroutine cross_prod( a, b, c )
    !> @brief Cross product of two vectors
    implicit none
    real(rp), dimension(3), intent(in) :: a
    real(rp), dimension(3), intent(in) :: b
    real(rp), dimension(3), intent(out) :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine cross_prod


  recursive function gcd_rec(u, v) result(gcd)
    !! greatest common denominator
    integer             :: gcd
    integer(ip), intent(in) :: u, v

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
    real(rp), dimension(3,3), intent(in) :: a, b
    real(rp), intent(out) :: dist

    integer(ip) :: i, j

    dist = 0.0_rp
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
    real(rp), dimension(3,3), intent(in) :: a, b
    real(rp), intent(out) :: dist

    real(rp), dimension(3,3) :: rij
    real(rp) :: tr

    rij = matmul( a, transpose(b) )
    rij = transpose(rij)
    tr = rij(1,1) + rij(2,2) + rij(3,3) - 3.0_rp
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
    interface
       subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
         use, intrinsic :: iso_fortran_env, only: dp=>real64
         integer, intent(out) :: info
         integer, intent(in) :: lda, ldvl, ldvr, lwork, n
         character(1), intent(in) :: jobvl, jobvr
         real(dp), intent(inout) :: a(lda, *), vl(ldvl, *), vr(ldvr, *), wi(*), wr(*)
         real(dp), intent(out) :: work(max(1,lwork))
         intrinsic :: max
       end subroutine dgeev
    end interface
    INTEGER(IP),              intent(in) :: n
    REAL(RP), DIMENSION(n,n), intent(inout) :: A
    REAL(RP), DIMENSION(n),   intent(out) :: eigvals
    INTEGER(IP),              intent(in) :: vec
    REAL(RP), DIMENSION(n) :: eigvals_i !! imaginary part of the eigenvalues
    REAL(RP), DIMENSION(n,n) :: eigvec
    INTEGER(IP) :: lda
    INTEGER(IP) :: lwork
    REAL(RP) :: Dummy(1000)
    INTEGER(IP) :: info
    CHARACTER(len=1) :: getvec

    real(rp), dimension(3,3) :: mat
    getvec = 'N'
    if( vec == 1 ) getvec='V'
    lda = n
    eigvals_i(:) = 0.0_rp
    eigvec(:,:) = 0.0_rp

    mat = A
    !! test workspace
    lwork = -1
    call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, dummy, 1, eigvec, n, dummy, lwork, info)

    !! choose optimal size of workspace (as in example from intel website)
    lwork = min( 1000, nint(dummy(1)) )
    !! compute stuffs
    call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, dummy, 1, eigvec, n, dummy, lwork, info)

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
    integer(ip) :: err
    character(len=10), intent(in) :: pg
    integer(ip), intent(in) :: Nop
    integer(ip) :: exp_Nop

    err = 0

    !! compute expected Nop based on pg string
    exp_Nop = get_expected_Nop( pg )
    if( exp_Nop < 1 ) then
       write(*,*) "at:",__FILE__, " line:",__LINE__
       return
    end if

    if( Nop < exp_Nop ) then
       !! if Nop is less than should be, return err=-1
       err=-1
    elseif( Nop > exp_Nop ) then
       !! if Nop is more than should be, return err=1
       err = 1
    end if

  end function check_pg_Nop


  function get_expected_Nop( pg ) result( Nop )
    !! compute expected number of symmetry operations, given a pg name.
    !! formula: Nop = order * multiplier * mod
    !! where:
    !! `order` is the order of the C,D,S groups, or: `order(T)=12`, `order(O)=24`, `order(I)=60`
    !! `multiplier=2` for D groups, `multiplier=1` otherwise
    !! `mod=2` if group name has any of: v, h, d, s, i as the last letter, and `mod=1` otherwise
    !!
    !! On error, return Nop = -1
    implicit none
    character(len=*), intent(in) :: pg
    integer(ip) :: Nop
    !!local
    character(len=1) :: letter, mod
    character(len=10) :: pg_cpy
    integer(ip) :: order, ntot
    integer(ip) :: multiply_mod, multiply_letter

    Nop = -1

    order = 1
    mod = "x"

    pg_cpy = pg
    ntot = len_trim(pg_cpy)
    !! empty or overflow string?
    if( ntot < 1 .or. ntot > 10 ) then
       write(*,*) "at:",__FILE__," line:",__LINE__
       write(*,*) "ERROR:: string `pg` has invalid length:",ntot
       write(*,*) "`pg` string is:", pg
       return
    end if

    !! read first letter
    read( pg_cpy, "(a1)") letter
    ntot = ntot - 1

    if( ntot > 0 ) then
       !! cut first letter from string
       pg_cpy = trim(pg(2:))
       !! check last letter for mod, if it is one of those, copy it
       select case( pg_cpy(ntot:ntot) )
       case( "v", "h", "d", "s", "i" )
          mod = pg_cpy(ntot:ntot)
          !! shift n back by 1 letter
          ntot = ntot - 1
       end select
    end if

    if( ntot > 0 ) then
       !! read order integer
       select case( letter )
       case( "C", "D", "S" )
          !! read 1:n for order
          read( pg_cpy(1:ntot), * ) order
       end select
    end if

    !! compute expected number of elements:
    multiply_letter = 1
    multiply_mod    = 1
    if( letter == "T" ) order = 12
    if( letter == "O" ) order = 24
    if( letter == "I" ) order = 60
    !! D groups have order*2 elements
    if( letter == "D" ) multiply_letter = 2
    !! any nonzero mod doubles the number of elements
    if( mod /= "x" ) multiply_mod = 2
    !! final number of operations
    Nop = order * multiply_letter * multiply_mod
  end function get_expected_Nop


  function find_inversion( nbas, op ) result( has_inversion )
    !! find if there is inversion in the list
    implicit none
    logical :: has_inversion
    integer(ip), intent(in) :: nbas
    character(len=1), dimension(nbas), intent(in) :: op
    integer(ip) :: i

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
    integer(ip), intent(in) :: nbas
    character(len=1), dimension(nbas), intent(in) :: op
    integer(ip), dimension(nbas), intent(in) :: n_int
    integer(ip) :: i

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
    integer(ip), intent(in) :: nbas
    character(len=1), dimension(nbas), intent(in) :: op
    integer(ip), dimension(nbas), intent(in) :: n_int
    integer(ip) :: i

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
    real(rp), dimension(3,3), intent(in) :: rmat
    real(rp), dimension(3), intent(in) :: bfield

    real(rp), dimension(3) :: bdir, vec
    real(rp) :: det, dotp

    is_valid = .false.

    !! direction of B field
    bdir = bfield/norm2(bfield)

    call determinant3x3(rmat, det)

    !! det(M)*MB
    vec = matmul( rmat, bdir )
    vec = vec * det

    !! compare vec to bdir, shoudl be equal: dot = 1.0
    dotp = dot_product( vec, bdir )
    if( dotp .gt. 0.999_rp ) is_valid = .true.

  end function op_valid_ext_b

  subroutine construct_reflection( ax, rmat )
    !! construct reflection matrix from axis, as:
    !! rmat = I - 2.0* ( <x x^T> / <x^T x>), where x is the axis
    real(rp), dimension(3), intent(in) :: ax
    real(rp), dimension(3,3), intent(out) :: rmat

    real(rp), dimension(3,3) :: id, rr
    integer(ip) :: i, j

    id(:,:) = 0.0_rp
    do i = 1, 3
       id(i,i) = 1.0_rp
    end do

    !! outer product
    forall (i=1:3)
       forall(j=1:3) rr(i,j) = ax(i)*ax(j)
    end forall
    rr(:,:) = rr / norm2(ax)

    rmat(:,:) = id(:,:) - 2.0_rp*rr(:,:)

  end subroutine construct_reflection
  subroutine construct_rotation(ax_in, angle, rmat)
    !! wikipedia/rotation_matrix#quaternion
    implicit none
    real(rp), dimension(3), intent(in) :: ax_in
    real(rp), intent(in) :: angle
    real(rp), dimension(3,3), intent(out) :: rmat

    real(rp), dimension(3) :: ax
    real(rp) :: c, s, c_one

    ax = ax_in/norm2(ax_in)

    c = cos(angle)
    s = sin(angle)
    c_one = 1.0_rp - c

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
    real(rp), intent(inout) :: ax(3)

    real(rp) :: flip

    flip = 1.0_rp
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
