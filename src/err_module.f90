module err_module

  implicit none
  private
  public :: get_err_msg

  integer, parameter, public :: &
       !!
       !! IRA errors
       ERR_TOO_SMALL_KMAX = -1, &
       ERR_OTHER          = -2, &
       !!
       !! SOFI errors
       ERR_DETERMINANT    = -3, &
       ERR_ACOS_ARG       = -4

  character(len=*), parameter :: warn = "WARNING FROM IRA/SOFI library :::"
  character(len=*), parameter :: err = "ERROR FROM IRA/SOFI library :::"
contains

  function get_err_msg( ierr )result(msg)
    integer, intent(in) :: ierr
    character(:), allocatable :: msg

    character(len=512) :: str


    select case( ierr )
    case( ERR_TOO_SMALL_KMAX )
       write(str, "(5x, a,1x,a)") err, &
            "No basis could be found. Possible cause: too small `kmax_factor` value."

    case( ERR_OTHER )
       write(str, "(5x,a,1x,a,1x,i0)") err, &
            "Some other error has occured ... ierr =",ierr

    case( ERR_DETERMINANT )
       write(str, "(5x,a,1x,a)") err, &
            "Value of matrix determinant out of range!"

    case( ERR_ACOS_ARG )
       write(str, "(5x,a,1x,a)") err, &
            "Invalid value for acos argument, should be strictly on range [-1:1]"

    case default
       write(str, "(5x, a,1x,a,1x,i0)") warn, &
            "Unknown error code:", ierr

    end select

    allocate( msg, source=trim(str) )
  end function get_err_msg


end module err_module

