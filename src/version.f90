
!> @details
!! get the current IRA version, string and date of release
!!
!! @param[out] string :: version string
!! @param[out] date :: date of release, format: YYYYmmdd
subroutine ira_get_version( string, date )
  implicit none
  character(len=5), intent(out) :: string
  integer, intent(out) :: date

  !! version string
  string = "2.2.0"


  !! date string, format: YYYYmmdd
  date = 20250911


end subroutine ira_get_version
