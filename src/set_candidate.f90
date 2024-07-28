!> @brief
!! Subroutine to set candidate central atoms for structures 1 and 2.
!!
!! The decision is currently based simply on the total number of atoms,
!! if the number of atoms is equal, then both candidates get a -1 value,
!! otherwise the value of candidates is the atomic index of desired
!! central atom.
!!
!! UPDATE NEEDED:
!! candidate central atoms in 2 can only be of the same typ as in 1
subroutine set_candidates( nat1, typ1, coords1, &
                           nat2, typ2, coords2, &
                           candidate1, candidate2 )

  use ira_precision
  implicit none
  integer(ip),                  intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1
  real(rp), dimension(3,nat1),  intent(in) :: coords1
  integer(ip),                  intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2
  real(rp), dimension(3,nat2),  intent(in) :: coords2
  integer(ip), dimension(nat1), intent(out) :: candidate1
  integer(ip), dimension(nat2), intent(out) :: candidate2

  integer(ip) :: i, dnat, k

  candidate1(:) = 0
  candidate2(:) = 0
  !!
  !! some preprocessing
  dnat = nat2-nat1
  !!
  if( dnat .eq. 0 ) then
     !! nat1 = nat2

     !! value -1 indicates some special vector (geometrical center)
     candidate1(1) = -1

     !! value -1 indicates some special vector (geometrical center)
     candidate2(1) = -1


  elseif( dnat .gt. 0) then
     !! nat2 > nat1

     !! in struc 1 take the first atom
     candidate1(1) = 1

     !! in struc 2 take all atoms
     k = 1
     do i = 1, nat2
        if( typ2(i) .ne. typ1(1) ) cycle
        candidate2(k) = i
        k = k + 1
     end do

  else
     write(*,*) 'error in set_candidate'
     stop
  endif



end subroutine set_candidates



subroutine select_rc( nat, coords, c_idx, rc )
  use ira_precision
  implicit none
  integer(ip), intent(in) :: nat
  real(rp), dimension(3,nat), intent(in) :: coords
  integer(ip), intent(in) :: c_idx
  real(rp), dimension(3), intent(out) :: rc

  if( c_idx .eq. 0 ) then
     write(*,*) 'ERROR in select_rc'
     return
  endif

  if( c_idx .eq. -1 ) then
     !! rc is geometric center
     rc = sum(coords(:,:),2)/nat
  else
     !! rc is vector of atom c_idx
     rc = coords(:,c_idx)
  endif

  return
end subroutine select_rc

