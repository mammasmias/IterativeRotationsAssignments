



module timer

  !> @details
  !! a simple timer module with 4 memory slots
  !! Call as:
  !!
  !! @code{.f90}
  !!   use timer
  !!   ! start timing slot 1
  !!   call timer_start( LOC_T1 )
  !!   ...
  !!   ! start timing slot 2
  !!   call timer_start( LOC_T2 )
  !!   ...
  !!   ! stop timer slot 2
  !!   call timer_stop( LOC_T2 )
  !!   ...
  !!   ...
  !!   ! stop timing slot 1
  !!   call timer_stop( LOC_T1 )
  !!
  !!   ! print timings
  !!   call timer_print()
  !! @endcode

  implicit none
  public


  !! slots for accessing the timer
  integer, parameter :: &
       LOC_T1 = 1, &
       LOC_T2 = 2, &
       LOC_T3 = 3, &
       LOC_T4 = 4


  !! the timer, store value dt1 which is clock when timer started
  type tim
     real :: dt1
   contains
     procedure :: start => tim_start
     procedure :: end => tim_end
  end type tim

  interface tim
     module procedure :: tim_constructor
  end interface tim


  type( tim ), pointer :: &
       t_1 => null(), &
       t_2 => null(), &
       t_3 => null(), &
       t_4 => null()
  real, protected :: timed(4)   !! array for total sum of each slot

contains

  function tim_constructor( ) result( this )
    type( tim ), pointer :: this
    allocate( tim::this )
  end function tim_constructor

  subroutine tim_start( self )
    implicit none
    class(tim), intent(inout) :: self
    !! save current clock into dt1
    call clock( self% dt1 )
  end subroutine tim_start

  function tim_end( self )result(dt)
    implicit none
    class(tim), intent(inout) :: self
    real :: dt
    real :: dt2
    !! get current clock
    call clock( dt2 )
    !! return difference of current clock and when i started the clock
    dt = dt2 - self% dt1
  end function tim_end



  subroutine timer_start( loc )
    !! start the timer for slot `loc`
    implicit none
    integer, intent(in) :: loc
    select case( loc )
    case( LOC_T1 ); if( .not. associated(t_1) ) t_1 => tim(); call t_1% start()
    case( LOC_T2 ); if( .not. associated(t_2) ) t_2 => tim(); call t_2% start()
    case( LOC_T3 ); if( .not. associated(t_3) ) t_3 => tim(); call t_3% start()
    case( LOC_T4 ); if( .not. associated(t_4) ) t_4 => tim(); call t_4% start()
    end select
  end subroutine timer_start

  subroutine timer_stop( loc )
    !! stop timer for slot `loc` and add dt to total sum of time for slot
    implicit none
    integer, intent(in) :: loc
    real :: dt
    select case( loc )
    case( LOC_T1 ); dt = t_1% end()
    case( LOC_T2 ); dt = t_2% end()
    case( LOC_T3 ); dt = t_3% end()
    case( LOC_T4 ); dt = t_4% end()
    end select
    timed( loc ) = timed(loc) + dt
  end subroutine timer_stop

  subroutine timer_print()
    !! print current total times for all slots
    implicit none
    write(*,"(1x,a)")       ":::: timer ::::"
    write(*,"(3x,a,2x,g0.8)") "t1  ::", timed(LOC_T1)
    write(*,"(3x,a,2x,g0.8)") "t2  ::", timed(LOC_T2)
    write(*,"(3x,a,2x,g0.8)") "t3  ::", timed(LOC_T3)
    write(*,"(3x,a,2x,g0.8)") "t4  ::", timed(LOC_T4)
    write(*,"(1x,a)")       ":::::::::::::::"
  end subroutine timer_print


  subroutine clock( time )
    use, intrinsic :: iso_fortran_env, only: int64
    real, intent(out) :: time

    integer( int64 ) :: tick, rate

    call system_clock( tick, rate )
    time = real( tick, kind(time) ) / real( rate, kind(time) )
  end subroutine clock

end module timer

