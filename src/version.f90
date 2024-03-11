
!> @details
!! get the current IRA version, string and date of release
!!
!! @param[out] string :: version string
!! @param[out] date :: date of release, format: YYYYmm
subroutine get_version( string, date )
  implicit none
  character(len=9), intent(out) :: string
  integer, intent(out) :: date

  !! version string
  string = "IRAv1.6.0"


  !! date string, format: YYYYmm
  date = 202403

end subroutine get_version
