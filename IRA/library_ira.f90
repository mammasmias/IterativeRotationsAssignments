
subroutine lib_cshda( nat1, typ1, coords1, nat2, typ2, coords2, thr, found, dists )bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat1
  type( c_ptr ), value, intent(in) :: typ1
  type( c_ptr ), value, intent(in) :: coords1
  integer(c_int), value, intent(in) :: nat2
  type( c_ptr ), value, intent(in) :: typ2
  type( c_ptr ), value, intent(in) :: coords2
  real( c_double ), value, intent(in) :: thr
  !!
  type( c_ptr ), intent(in) :: found
  type( c_ptr ), intent(in) :: dists

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_found
  type( c_ptr ), dimension(:), pointer :: p_atm1, p_atm2
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:), pointer :: p_dists

  !! local f memory
  real(c_double), dimension(3,3) :: f_lat
  integer :: i
  real(c_double) :: some_thr


  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_atm1, [nat1] )
  call c_f_pointer( coords2, p_atm2, [nat2] )
  call c_f_pointer( p_atm1(1), p_coords1, [3,nat1] )
  call c_f_pointer( p_atm2(1), p_coords2, [3,nat2] )

  call c_f_pointer( found, p_found, [nat2] )
  call c_f_pointer( dists, p_dists, [nat2] )

  some_thr=thr

  call cshda( nat1, p_typ1, p_coords1, nat2, p_typ2, p_coords2, some_thr, &
       p_found, p_dists )

end subroutine lib_cshda



subroutine lib_cshda_pbc( nat1, typ1, coords1, nat2, typ2, coords2, lat, thr, found, dists )bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), value, intent(in) :: nat1
  type( c_ptr ), value, intent(in) :: typ1
  type( c_ptr ), value, intent(in) :: coords1
  integer(c_int), value, intent(in) :: nat2
  type( c_ptr ), value, intent(in) :: typ2
  type( c_ptr ), value, intent(in) :: coords2
  type( c_ptr ), value, intent(in) :: lat
  real( c_double ), value, intent(in) :: thr
  !!
  type( c_ptr ), intent(in) :: found
  type( c_ptr ), intent(in) :: dists

  !! f ptrs
  integer(c_int), dimension(:), pointer :: p_typ1, p_typ2, p_found
  type( c_ptr ), dimension(:), pointer :: p_atm1, p_atm2
  real( c_double ), dimension(:,:), pointer :: p_coords1, p_coords2
  real( c_double ), dimension(:), pointer :: p_dists, p_lat

  !! local f memory
  real(c_double), dimension(3,3) :: f_lat
  integer :: i
  real(c_double) :: some_thr


  !! connect c ptrs to f
  call c_f_pointer( typ1, p_typ1, [nat1] )
  call c_f_pointer( typ2, p_typ2, [nat2] )
  call c_f_pointer( coords1, p_atm1, [nat1] )
  call c_f_pointer( coords2, p_atm2, [nat2] )
  call c_f_pointer( p_atm1(1), p_coords1, [3,nat1] )
  call c_f_pointer( p_atm2(1), p_coords2, [3,nat2] )
  call c_f_pointer( lat, p_lat, [9] )
  f_lat = reshape( p_lat, [3,3] )
  f_lat = transpose( f_lat )


  call c_f_pointer( found, p_found, [nat2] )
  call c_f_pointer( dists, p_dists, [nat2] )

  some_thr=thr

  call cshda_pbc( nat1, p_typ1, p_coords1, nat2, p_typ2, p_coords2, f_lat, &
       some_thr, p_found, p_dists )

end subroutine lib_cshda_pbc

