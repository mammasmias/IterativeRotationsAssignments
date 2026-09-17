#define ASSERT( cond, msg )\
  if( .not. cond ) then; \
     write(*,"('> error at >:',1x,a,'::',i0,':')")__FILE__,__LINE__ ;\
     write(*,"('>',1x,*(g0))") msg ;\
     error stop 1 ;\
  end if

module dbg
  use ira_precision, only: ira_rp => rp
  private
  public :: tostr
  public :: within
  public :: check_permutation
  public :: random_permutation
  public :: inverse_perm

  interface tostr
     procedure :: tostr_r0, tostr_r1
  end interface tostr

contains

  function tostr_r0( val, fmt )result(str)
    class(*), intent(in) :: val
    character(*), intent(in),optional :: fmt
    character(:), allocatable :: fmt1
    character(:), allocatable :: str
    character(:), allocatable :: long_str
    fmt1="(g0)"
    if(present(fmt))fmt1=fmt
    select type(val)
    type is(integer)
       allocate(character(len=2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(real)
       allocate(character(len=2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(real(ira_rp))
       allocate(character(len=2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(logical)
       allocate(character(len=2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(character(*))
       allocate(character(len=32*len(val)) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    end select
    str=trim(long_str)
  end function tostr_r0
  function tostr_r1( val, fmt )result(str)
    class(*), intent(in) :: val(:)
    character(*), intent(in), optional :: fmt
    character(:), allocatable :: fmt1
    character(:), allocatable :: str
    character(:), allocatable :: long_str
    fmt1="(*(g0,:,1x))"
    if(present(fmt))fmt1=fmt
    select type(val)
    type is(integer)
       allocate(character(len=size(val)*2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(real)
       allocate(character(len=size(val)*2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(real(ira_rp))
       allocate(character(len=size(val)*2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is(logical)
       allocate(character(len=size(val)*2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    type is( character(*) )
       allocate(character(len=size(val)*2048) :: long_str)
       write(long_str, "("//trim(fmt1)//")" )val
    end select
    str=trim(long_str)
  end function tostr_r1



  elemental function within( a, b, tol )result(res)
    class(*), intent(in) :: a, b
    real, intent(in) :: tol
    logical :: res
    select type(a)
    type is( integer     ); select type(b); type is(integer); res = abs( a - b ) <= tol; end select
    type is( real        ); select type(b); type is(real) ; res = abs( a - b ) <= tol; end select
    type is( real(ira_rp)); select type(b); type is(real(ira_rp)) ; res = abs( a - b ) <= tol; end select
    end select
  end function within


  function check_permutation( permutation )result(ierr)
    implicit none
    integer, intent(in) :: permutation(:)
    integer :: ierr
    integer :: perm_seen(size(permutation))
    integer :: n, ip
    ierr = -1
    n = size(permutation)
    perm_seen = 0

    do ip = 1, n
       if( permutation(ip) < 1 .or. permutation(ip) > n ) then
          write(*,*) "ERROR in check_permutation"
          write(*,"(1x,a,i0)") "permutation entry outside the range 1:",n
          write(*,"(1x,a,i0)") "permutation(i) = ",permutation(ip)
          return
       end if
       perm_seen( permutation(ip) ) = perm_seen( permutation(ip) ) + 1
    end do
    if( any( perm_seen /= 1 ) ) then
       write(*,*) "ERROR in check_permutation"
       write(*,*) "permutation is not a bijection:"
       write(*,"(10(i3,1x))") permutation
       return
    end if
    ierr = 0
  end function check_permutation

  subroutine random_permutation( n, list )
    !! generate a list of random indices of size n, such that
    !! no index repeats
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(out) :: list

    integer :: i, idx
    real :: z
    logical :: old

    !! initial values
    list(:) = 0

    do i = 1, n
       old = .true.
       do while( old )
          call random_number(z)
          !! generate index randomly in the range [1:n]
          idx = int( z*n ) + 1
          !! if this index already in list, skip and geenrate new
          if( any(list .eq. idx ) ) cycle
          !! if not, add it to list and stop loop for current index
          list(i) = idx
          old = .false.
       end do
    end do
  end subroutine random_permutation

  function inverse_perm( p )result(ip)
    !! The inverse of a permutation is equal to index of each value,
    !! e.g. ip(1) is the index of value 1 in p array.
    !! `p` on input needs to be a valid permutation (no invalid values, or zeros!)
    implicit none
    integer, intent(in) :: p(:)
    integer :: ip(size(p))
    integer :: i
    do i = 1, size(p)
       ip(i) = findloc(p,i,dim=1)
    end do
  end function inverse_perm

end module dbg

