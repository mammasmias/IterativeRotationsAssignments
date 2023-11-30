module err_module

  implicit none
  private
  public :: get_err_msg

  integer, parameter, public :: &
       ERR_TOO_SMALL_KMAX = -1, &
       ERR_OTHER          = -2

  character(len=*), parameter :: warn = "WARNING FROM IRA/SOFI library :::"
  character(len=*), parameter :: err = "ERROR FROM IRA/SOFI library :::"
contains

  function get_err_msg( ierr )result(msg)
    integer, intent(in) :: ierr
    character(:), allocatable :: msg

    character(len=512) :: line

    select case( ierr )
    case( ERR_TOO_SMALL_KMAX )
       write(line, "(5x, a,1x,a)") err,&
            "No basis could be found. Possible cause: too small `kmax_factor` value."

    case( ERR_OTHER )
       write(line, "(5x,a,1x,a,1x,i0)") err, &
            "Some other error has occured ... ierr =",ierr

    case default
       write(line, "(5x, a,1x,a,1x,i0)") warn, &
            "Unknown error code:", ierr

    end select

    allocate( msg, source=trim(line) )
  end function get_err_msg


end module err_module
