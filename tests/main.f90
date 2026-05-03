#include "dbg.f90"
program main
  use dbg
  use ira_precision, only: ira_rp => rp
  implicit none
  real(ira_rp) :: r
  character(len=5) :: str
  integer :: int
  r = 1.0_ira_rp
  ASSERT( within(r, 1.0_ira_rp, 1e-2), "r is: "//tostr(r) )
  call ira_get_version( str, int )

  block
    ! test SOFI, this calls lapack for diagonalization in analmat
    real(ira_rp) :: rmat(3,3), rmat_expected(3,3), ax_in(3), angle_in
    integer :: ierr
    character(len=1) :: op
    integer :: n, p
    real(ira_rp) :: ax(3), angle, ax_expected(3)
    ! specify a rotation operation from axis and angle
    angle_in = real(1.0/3.0, kind=ira_rp)
    ax_in = real([1.0, 1.0, 1.0], kind=ira_rp)
    call sofi_construct_operation( "C", ax_in, angle_in, rmat, ierr )
    ASSERT( ierr == 0, "Expected ierr=0, got ierr="//tostr(ierr) )
    rmat_expected(1,:) = [ 0.0_ira_rp, 0.0_ira_rp, 1.0_ira_rp ]
    rmat_expected(2,:) = [ 1.0_ira_rp, 0.0_ira_rp, 0.0_ira_rp ]
    rmat_expected(3,:) = [ 0.0_ira_rp, 1.0_ira_rp, 0.0_ira_rp ]
    ASSERT( all(within(rmat,rmat_expected, 1e-2)), "rmat is off?" )

    ! analyse the matrix to obtain axis, angle, etc
    call sofi_analmat( rmat, op, n, p, ax, angle, ierr )
    ASSERT( ierr == 0, "Expected ierr=0, got ierr="//tostr(ierr) )
    ASSERT( op=="C", "Expected op='C', got op="//op )
    ASSERT( n == 3, "Expected n=3, got n="//tostr(n) )
    ASSERT( p == 1, "Expected p=1, got p="//tostr(p) )
    ax_expected = ax_in/sqrt(3.0_ira_rp)
    ASSERT( all(within(ax, ax_expected, 1e-6)), "Expected ax="//tostr(ax_expected,fmt="(3(g0.4,:,1x))")//", got ax="//tostr(ax,fmt="(3(g0.4,:,1x))"))
    ASSERT( within(angle, angle_in, 1e-12), "Expected angle="//tostr(angle_in)//", got angle="//tostr(angle))
  end block

end program
