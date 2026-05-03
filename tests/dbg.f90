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

end module dbg

