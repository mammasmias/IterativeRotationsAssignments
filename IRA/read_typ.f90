!!
!! Copyright (C) 2021, MAMMASMIAS Consortium
!! Written by: Miha Gunde
!!
!! This file is distributed under the terms of GNU General Public License.
!! See the license at: http://www.gnu.org/licenses/gpl-3.0.txt
!!
module read_typ

  !! total number of atomic types
  integer :: ntyp = 0

  !! array of all known atomic types, in string
  character(len=2), allocatable :: all_ctyp(:)

contains

    integer function get_int_typ( this_ctyp ) result( ityp )
      !! this function transforms the atomic type from character to an integer,
      !! according to array of all known types all_ctyp(:)
      character(len=2), intent(in) :: this_ctyp

      integer :: k

      !! if this function is called first time, allocate 100 typs
      if( .not. allocated(all_ctyp) ) then
         allocate( all_ctyp(1:100),source='')
      endif

      if( .not.any(all_ctyp(:) .eq. this_ctyp) ) then
         !! is new typ
         ntyp = ntyp + 1
         all_ctyp( ntyp ) = this_ctyp
         ityp = ntyp
         return
      else
         !! is known typ
         do k = 1, ntyp
            if( all_ctyp(k) == this_ctyp ) exit
         end do
         ityp = k
         return
      endif
    end function get_int_typ


end module read_typ
