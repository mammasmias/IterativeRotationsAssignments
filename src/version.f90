#ifdef IRA_STRVERSION
#undef IRA_STRVERSION
#endif
#ifdef IRA_DATEVERSION
#undef IRA_DATEVERSION
#endif

#define IRA_STRVERSION STRVERSION
#define IRA_DATEVERSION DATEVERSION

!> @details
!! get the current IRA version, string and date of release
!!
!! @param[out] string :: version string
!! @param[out] date :: date of release, format: YYYYmmdd
subroutine ira_get_version( string, date )
  implicit none
  character(len=5), intent(out) :: string
  integer, intent(out) :: date

  character(*), parameter :: strversion=IRA_STRVERSION
  integer, parameter :: dateversion=IRA_DATEVERSION
  !! version string
  string = strversion


  !! date string, format: YYYYmmdd
  date = dateversion


end subroutine ira_get_version
