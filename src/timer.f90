



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

  use ira_precision
  implicit none
  public


  !! slots for accessing the timer
  integer(ip), parameter :: &
       LOC_T1 = 1, &
       LOC_T2 = 2, &
       LOC_T3 = 3, &
       LOC_T4 = 4


  !! a single timer, store value dt1 which is clock when timer started
  type tim
     real(rp) :: dt1
   contains
     procedure :: start => tim_start
     procedure :: end => tim_end
  end type tim

  interface tim
     module procedure :: tim_constructor
  end interface tim


  !! definition of local timer for use within single routine
  type local_timer
     type( tim ), pointer :: &
          t_1 => null(), &
          t_2 => null(), &
          t_3 => null(), &
          t_4 => null()
     real(rp) :: timed(4)   !! array for total sum of each slot
     character(:), allocatable :: tag1, tag2, tag3, tag4
   contains
     procedure :: start => timer_start
     procedure :: stop => timer_stop
     procedure :: tag => timer_tag
     procedure :: print => timer_print
  end type local_timer


  !! definiton of global_timer which works accross routines ... not really
  type( tim ), pointer :: &
       t_1 => null(), &
       t_2 => null(), &
       t_3 => null(), &
       t_4 => null()
  real(rp), protected :: timed(4)   !! array for total sum of each slot
  character(:), allocatable :: tag1, tag2, tag3, tag4


contains

  !!===== functions for single value of timer ===
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
    real(rp) :: dt
    real(rp) :: dt2
    !! get current clock
    call clock( dt2 )
    !! return difference of current clock and when i started the clock
    dt = dt2 - self% dt1
  end function tim_end
  !! ============



  !! ==== functions for local timer =======
  subroutine timer_start( self, loc )
    !! start the timer for slot `loc`
    implicit none
    class( local_timer ), intent(inout) :: self
    integer(ip), intent(in) :: loc
    select case( loc )
    case( LOC_T1 ); if( .not. associated(self% t_1) ) self% t_1 => tim(); call self% t_1% start()
    case( LOC_T2 ); if( .not. associated(self% t_2) ) self% t_2 => tim(); call self% t_2% start()
    case( LOC_T3 ); if( .not. associated(self% t_3) ) self% t_3 => tim(); call self% t_3% start()
    case( LOC_T4 ); if( .not. associated(self% t_4) ) self% t_4 => tim(); call self% t_4% start()
    end select
  end subroutine timer_start

  subroutine timer_stop( self, loc )
    !! stop timer for slot `loc` and add dt to total sum of time for slot
    implicit none
    class( local_timer ), intent(inout) :: self
    integer(ip), intent(in) :: loc
    real(rp) :: dt
    select case( loc )
    case( LOC_T1 ); dt = self% t_1% end()
    case( LOC_T2 ); dt = self% t_2% end()
    case( LOC_T3 ); dt = self% t_3% end()
    case( LOC_T4 ); dt = self% t_4% end()
    end select
    self% timed( loc ) = self% timed(loc) + dt
  end subroutine timer_stop

  subroutine timer_tag( self, loc, tag )
    !! set a tag to loc
    implicit none
    class( local_timer ), intent(inout) :: self
    integer(ip), intent(in) :: loc
    character(*), intent(in) :: tag
    select case( loc )
    case( LOC_T1 ); self% tag1 = tag
    case( LOC_T2 ); self% tag2 = tag
    case( LOC_T3 ); self% tag3 = tag
    case( LOC_T4 ); self% tag4 = tag
    end select
  end subroutine timer_tag

  subroutine timer_print( self )
    !! print current total times for all slots
    implicit none
    class( local_timer ), intent(inout) :: self
    write(*,"(1x,a)")       ":::: local timer ::::"
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t1  ::", self% tag1, self% timed(LOC_T1)
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t2  ::", self% tag2, self% timed(LOC_T2)
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t3  ::", self% tag3, self% timed(LOC_T3)
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t4  ::", self% tag4, self% timed(LOC_T4)
    write(*,"(1x,a)")       ":::::::::::::::"
  end subroutine timer_print
  !!===========================


  !! ========= functions for global timer ==============
  subroutine global_timer_start( loc )
    !! start the timer for slot `loc`
    implicit none
    integer(ip), intent(in) :: loc
    select case( loc )
    case( LOC_T1 ); if( .not. associated(t_1) ) t_1 => tim(); call t_1% start()
    case( LOC_T2 ); if( .not. associated(t_2) ) t_2 => tim(); call t_2% start()
    case( LOC_T3 ); if( .not. associated(t_3) ) t_3 => tim(); call t_3% start()
    case( LOC_T4 ); if( .not. associated(t_4) ) t_4 => tim(); call t_4% start()
    end select
  end subroutine global_timer_start

  subroutine global_timer_stop( loc )
    !! stop timer for slot `loc` and add dt to total sum of time for slot
    implicit none
    integer(ip), intent(in) :: loc
    real(rp) :: dt
    select case( loc )
    case( LOC_T1 ); dt = t_1% end()
    case( LOC_T2 ); dt = t_2% end()
    case( LOC_T3 ); dt = t_3% end()
    case( LOC_T4 ); dt = t_4% end()
    end select
    timed( loc ) = timed(loc) + dt
  end subroutine global_timer_stop

  subroutine global_timer_tag( loc, tag )
    !! set a tag to loc
    implicit none
    integer(ip), intent(in) :: loc
    character(*), intent(in) :: tag
    select case( loc )
    case( LOC_T1 ); tag1 = tag
    case( LOC_T2 ); tag2 = tag
    case( LOC_T3 ); tag3 = tag
    case( LOC_T4 ); tag4 = tag
    end select
  end subroutine global_timer_tag

  subroutine global_timer_print( )
    !! print current total times for all slots
    implicit none
    write(*,"(1x,a)")       ":::: global timer ::::"
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t1  ::", tag1, timed(LOC_T1)
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t2  ::", tag2, timed(LOC_T2)
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t3  ::", tag3, timed(LOC_T3)
    write(*,"(3x,a,2x, 'tag:',1x,a10, 2x, g0.8)") "t4  ::", tag4, timed(LOC_T4)
    write(*,"(1x,a)")       ":::::::::::::::"
  end subroutine global_timer_print
  !! =========================




  subroutine clock( time )
    use, intrinsic :: iso_fortran_env, only: int64
    real(rp), intent(out) :: time

    integer( int64 ) :: tick, rate

    call system_clock( tick, rate )
    time = real( tick, kind(time) ) / real( rate, kind(time) )
  end subroutine clock

end module timer

