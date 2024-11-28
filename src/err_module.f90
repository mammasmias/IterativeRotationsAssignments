module err_module

  use ira_precision
  implicit none
  private
  public :: get_err_msg

  integer(ip), parameter, public :: &
       !!
       !! IRA errors
       ERR_TOO_SMALL_KMAX = -1, &
       ERR_OTHER          = -2, &
       ERR_BETA           = -3, &
       ERR_SVD            = -9, &
       !!
       !! SOFI errors
       ERR_DETERMINANT    = -4, &
       ERR_ACOS_ARG       = -5, &
       ERR_SIZE_NMAX      = -6, &
       ERR_LIST_TOO_SMALL = -7, &
       ERR_UNKNOWN_OP     = -8

  character(len=*), parameter :: warn = "WARNING FROM IRA/SOFI library :::"
  character(len=*), parameter :: err = "ERROR FROM IRA/SOFI library :::"
contains

  function get_err_msg( ierr )result(msg)
    integer(ip), intent(in) :: ierr
    character(:), allocatable :: msg

    character(len=512) :: str


    select case( ierr )
    case( ERR_TOO_SMALL_KMAX )
       write(str, "(5x, a,1x,a)") err, &
            "No basis could be found. Possible cause: too small `kmax_factor` value."

    case( ERR_OTHER )
       write(str, "(5x,a,1x,a,1x,i0)") err, &
            "Some other error has occured ... ierr =",ierr

    case( ERR_BETA )
       write(str, "(5x,a,1x,a)") err, &
            "Cannot set beta basis!"

    case( ERR_SVD )
       write(str, "(5x,a,1x,a)") err, &
            "Lapack could not compute SVD."

    case( ERR_DETERMINANT )
       write(str, "(5x,a,1x,a)") err, &
            "Value of matrix determinant out of range!"

    case( ERR_ACOS_ARG )
       write(str, "(5x,a,1x,a)") err, &
            "Invalid value for acos argument, should be strictly on range [-1:1]"

    case( ERR_SIZE_NMAX )
       write(str, "(5x,a,1x,a)") err, &
            "The size of list for symmetry operation matrices is not nmax!"

    case( ERR_LIST_TOO_SMALL )
       write(str, "(5x,a,1x,a)") err, &
            "Number of symm operations exceeds the list size! Check your configuration."

    case( ERR_UNKNOWN_OP )
       write(str, "(5x,a,1x,a)") err, &
            "Unknown operation!"

    case default
       write(str, "(5x, a,1x,a,1x,i0)") warn, &
            "Unknown error code:", ierr

    end select

    allocate( msg, source=trim(str) )
  end function get_err_msg


end module err_module

