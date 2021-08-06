!!
!! Copyright (C) 2021, MAMMASMIAS Consortium
!! Written by: Miha Gunde
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!=====================================================
!! Description: 
!! This module reads atomic type and converts it to an appropriate integer.
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
