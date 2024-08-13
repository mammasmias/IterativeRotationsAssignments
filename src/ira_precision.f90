module ira_precision
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  private
  public :: rp, ip

  !! equivalent to c_int and c_double
  integer, parameter :: ip = int32
  integer, parameter :: rp = real64
end module ira_precision
