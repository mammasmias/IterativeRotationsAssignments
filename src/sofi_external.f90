

!! routines to be called as external (without loading ira_mod)

subroutine sofi_compute_all( nat, typ, coords, sym_thr, prescreen_ih, &
     nmat, mat_list, perm_list, &
     op_list, n_list, p_list, &
     ax_list, angle_list, dHausdorff_list, pg, n_prin_ax, prin_ax, &
     ierr )
  use ira_precision
  use m_sofi_tools, only: nmax
  use ira_mod, only: sofi_compute_all_x=>sofi_compute_all
  implicit none
  integer(ip),                 intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ
  real(rp), dimension(3,nat),  intent(in) :: coords
  real(rp),                    intent(in) :: sym_thr
  logical,                 intent(in) :: prescreen_ih
  integer(ip),                       intent(out) :: nmat
  real(rp), dimension(3,3,nmax),     intent(out) :: mat_list
  integer(ip), dimension(nat, nmax), intent(out) :: perm_list
  character(len=1), dimension(nmax), intent(out) :: op_list
  integer(ip), dimension(nmax),      intent(out) :: n_list
  integer(ip), dimension(nmax),      intent(out) :: p_list
  real(rp), dimension(3, nmax),      intent(out) :: ax_list
  real(rp), dimension(nmax),         intent(out) :: angle_list
  real(rp), dimension(nmax),         intent(out) :: dHausdorff_list
  character(len=10),             intent(out) :: pg
  integer(ip),                       intent(out) :: n_prin_ax
  real(rp), dimension(3,nmax),       intent(out) :: prin_ax
  integer(ip),                       intent(out) :: ierr
  call sofi_compute_all_x( nat, typ, coords, sym_thr, prescreen_ih, &
       nmat, mat_list, perm_list, &
       op_list, n_list, p_list, &
       ax_list, angle_list, dHausdorff_list, pg, n_prin_ax, prin_ax, &
       ierr )
end subroutine sofi_compute_all

subroutine sofi_struc_pg( nat, typ_in, coords_in, sym_thr, pg, verb )
  use ira_precision
  use m_sofi_routines, only: sofi_struc_pg_x=>sofi_struc_pg
  implicit none
  integer(ip),                 intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ_in
  real(rp), dimension(3,nat),  intent(in) :: coords_in
  real(rp),                    intent(in) :: sym_thr
  character(len=10),       intent(out) :: pg
  logical,                 intent(in) :: verb
  call sofi_struc_pg_x( nat, typ_in, coords_in, sym_thr, pg, verb )
end subroutine sofi_struc_pg

subroutine sofi_get_symmops( nat, typ_in, coords_in, sym_thr, prescreen_ih, n_so, op_list, ierr )
  use ira_precision
  use m_sofi_tools, only: nmax
  use m_sofi_routines, only: sofi_get_symmops_x=>sofi_get_symmops
  implicit none
  integer(ip),                 intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ_in
  real(rp), dimension(3,nat),  intent(in) :: coords_in
  real(rp),                    intent(in) :: sym_thr
  logical,                 intent(in) :: prescreen_ih
  integer(ip),                           intent(out) :: n_so
  real(rp), dimension(3,3,nmax),         intent(out) :: op_list
  integer(ip),                           intent(out) :: ierr
  call sofi_get_symmops_x( nat, typ_in, coords_in, sym_thr, prescreen_ih, n_so, op_list, ierr )
end subroutine sofi_get_symmops

subroutine sofi_get_perm( nat, typ, coords, nbas, bas_list, perm_list, dHausdorff_list )
  use ira_precision
  use m_sofi_routines, only: sofi_get_perm_x=>sofi_get_perm
  implicit none
  integer(ip),                   intent(in) :: nat
  integer(ip), dimension(nat),   intent(in) :: typ
  real(rp), dimension(3,nat),    intent(in) :: coords
  integer(ip),                   intent(in) :: nbas
  real(rp), dimension(3,3,nbas), intent(in) :: bas_list
  integer(ip), dimension(nat, nbas), intent(out) :: perm_list
  real(rp), dimension(nbas),     intent(out) :: dHausdorff_list
  call sofi_get_perm_x( nat, typ, coords, nbas, bas_list, perm_list, dHausdorff_list )
end subroutine sofi_get_perm

subroutine sofi_get_combos( nat, typ, coords, nbas, bas_list, ierr )
  use ira_precision
  use m_sofi_tools, only: nmax
  use m_sofi_routines, only: sofi_get_combos_x=>sofi_get_combos
  implicit none
  integer(ip), intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ
  real(rp), dimension(3, nat), intent(in) :: coords
  integer(ip), intent(inout) :: nbas
  real(rp), dimension(3, 3, nmax), intent(inout) :: bas_list
  integer(ip), intent(out) :: ierr

  call sofi_get_combos_x( nat, typ, coords, nbas, bas_list, ierr )
end subroutine sofi_get_combos

subroutine sofi_get_pg( nbas, op_list, pg, n_prin_ax, prin_ax, verb, ierr )
  use ira_precision
  use m_sofi_routines, only: sofi_get_pg_x => sofi_get_pg
  implicit none
  integer(ip),                   intent(in) :: nbas
  real(rp), dimension(3,3,nbas), intent(in) :: op_list
  character(len=10),         intent(out) :: pg
  integer(ip),                   intent(out) :: n_prin_ax
  real(rp), dimension(3,nbas),   intent(out) :: prin_ax
  logical,                   intent(in) :: verb
  integer(ip),                   intent(out) :: ierr
  call sofi_get_pg_x( nbas, op_list, pg, n_prin_ax, prin_ax, verb, ierr )
end subroutine sofi_get_pg

subroutine sofi_analmat( rmat, op, n, p, ax, angle, ierr )
  use ira_precision
  use m_sofi_routines, only: sofi_analmat_x=>sofi_analmat
  implicit none
  real(rp), dimension(3,3), intent(in) :: rmat
  character(len=1),     intent(out) :: op
  integer(ip),              intent(out) :: n
  integer(ip),              intent(out) :: p
  real(rp), dimension(3),   intent(out) :: ax
  real(rp),                 intent(out) :: angle
  integer(ip),              intent(out) :: ierr
  call sofi_analmat_x( rmat, op, n, p, ax, angle, ierr )
end subroutine sofi_analmat

subroutine sofi_ext_Bfield( n_op, op_list, b_field )
  use ira_precision
  use m_sofi_routines, only: sofi_ext_Bfield_x=>sofi_ext_Bfield
  implicit none
  integer(ip),                           intent(inout) :: n_op
  real(rp), dimension(1:3, 1:3, 1:n_op), intent(inout) :: op_list
  real(rp), dimension(3),                intent(in) :: b_field
  call sofi_ext_Bfield_x( n_op, op_list, b_field )
end subroutine sofi_ext_Bfield

subroutine sofi_construct_operation( op, axis, angle, matrix, ierr )
  use ira_precision
  use m_sofi_routines, only: sofi_construct_operation_x=>sofi_construct_operation
  implicit none
  character(len=1),     intent(in) :: op
  real(rp), dimension(3),   intent(in) :: axis
  real(rp),                 intent(in) :: angle
  real(rp), dimension(3,3), intent(out) :: matrix
  integer(ip),              intent(out) :: ierr
  call sofi_construct_operation_x( op, axis, angle, matrix, ierr )
end subroutine sofi_construct_operation

subroutine sofi_mat_combos( n_in, mat_in, n_out, mat_out )
  use ira_precision
  use m_sofi_routines, only: sofi_mat_combos_x=>sofi_mat_combos
  use m_sofi_tools, only: nmax
  implicit none
  integer(ip), intent(in) :: n_in
  real(rp), dimension(3,3,n_in), intent(in) :: mat_in
  integer(ip), intent(out) :: n_out
  real(rp), dimension(3,3,nmax), intent(out) :: mat_out
  call sofi_mat_combos_x( n_in, mat_in, n_out, mat_out )
end subroutine sofi_mat_combos

subroutine sofi_check_collinear( nat, coords, collinear, ax_o )
  use ira_precision
  use m_sofi_routines, only: sofi_check_collinear_x=>sofi_check_collinear
  implicit none
  integer(ip),                intent(in) :: nat
  real(rp), dimension(3,nat), intent(in) :: coords
  logical,                intent(out) :: collinear
  real(rp), dimension(3),     intent(out) :: ax_o
  call sofi_check_collinear_x( nat, coords, collinear, ax_o )
end subroutine sofi_check_collinear

subroutine sofi_get_err_msg( ierr, msg )
  use ira_precision
  use m_sofi_routines, only: sofi_get_err_msg_x => sofi_get_err_msg
  implicit none
  integer(ip), intent(in) :: ierr
  character(len=128), intent(out) :: msg
  call sofi_get_err_msg_x( ierr, msg )
end subroutine sofi_get_err_msg
