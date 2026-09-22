

!! routines to be called as external (without loading ira_mod)

subroutine ira_get_version(string, date)
  use ira_mod, only: ira_get_version_x => ira_get_version
  implicit none
  character(len=5), intent(out) :: string
  integer, intent(out) :: date
  call ira_get_version_x( string, date )
end subroutine ira_get_version

subroutine cshda( nat1, typ1, coords1, &
     nat2, typ2, coords2, &
     some_threshold, found, dists )
  use ira_precision
  use ira_mod, only: cshda_x => cshda
  implicit none
  integer(ip),                  intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1
  real(rp), dimension(3,nat1),  intent(in) :: coords1
  integer(ip),                  intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2
  real(rp), dimension(3,nat2),  intent(in) :: coords2
  real(rp),                     intent(in) :: some_threshold
  integer(ip), dimension(nat2), intent(out) :: found
  real(rp), dimension(nat2),    intent(out) :: dists
  call cshda_x( nat1, typ1, coords1, nat2, typ2, coords2, some_threshold, found, dists )
end subroutine cshda

subroutine cshda_from_cost( n2, n1, cost, found, dists )
  use ira_precision
  use ira_mod, only: cshda_from_cost_x => cshda_from_cost
  implicit none
  integer(ip), intent(in) :: n2
  integer(ip), intent(in) :: n1
  real(rp), intent(in) :: cost(n2,n1)
  integer(ip), intent(out) :: found(n2)
  real(rp), intent(out) :: dists(n2)
  call cshda_from_cost_x( n2, n1, cost, found, dists )
end subroutine cshda_from_cost

subroutine cshda_pbc( nat1, typ1, coords1, &
     nat2, typ2, coords2, lat2, &
     some_thr, found, dists )
  use ira_precision
  use ira_mod, only: cshda_pbc_x => cshda_pbc
  implicit none
  integer(ip),                  intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1
  real(rp), dimension(3,nat1),  intent(in) :: coords1
  integer(ip),                  intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2
  real(rp), dimension(3,nat2),  intent(in) :: coords2
  real(rp), dimension(3,3), intent(in) :: lat2
  real(rp),                     intent(in) :: some_thr
  integer(ip), dimension(nat2), intent(out) :: found
  real(rp), dimension(nat2),    intent(out) :: dists
  call cshda_pbc_x( nat1, typ1, coords2, nat2, typ2, coords2, lat2, some_thr, found, dists )
end subroutine cshda_pbc


subroutine ira_unify(nat1, typ1_in, coords1_in, candidate_1, &
     nat2, typ2_in, coords2_in, candidate_2, &
     kmax_factor, rotation, translation, permutation, hd_out, ierr)

  use ira_precision
  use ira_mod, only: ira_unify_x=> ira_unify
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1_in
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), dimension(nat1), intent(in) :: candidate_1
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2_in
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  integer(ip), dimension(nat2), intent(in) :: candidate_2
  real(rp), intent(in) :: kmax_factor
  real(rp), dimension(3, 3), intent(out) :: rotation
  real(rp), dimension(3), intent(out) :: translation
  integer(ip), dimension(nat2), intent(out) :: permutation
  real(rp), intent(out) :: hd_out
  integer(ip), intent(out) :: ierr
  call ira_unify_x( nat1, typ1_in, coords1_in, candidate_1, &
       nat2, typ2_in, coords2_in, candidate_2, &
       kmax_factor, rotation, translation, permutation, hd_out, ierr )
end subroutine ira_unify

subroutine ira_svd(nat1, typ1_in, coords1_in, &
     nat2, typ2_in, coords2_in, &
     kmax_factor, rotation, translation, permutation, &
     hd, rmsd, ierr)
  use ira_precision
  use ira_mod, only: ira_svd_x=>ira_svd
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1_in
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2_in
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  real(rp), intent(in) :: kmax_factor
  real(rp), dimension(3, 3), intent(out) :: rotation
  real(rp), dimension(3), intent(out) :: translation
  integer(ip), dimension(nat2), intent(out) :: permutation
  real(rp), intent(out) :: hd
  real(rp), intent(out) :: rmsd
  integer(ip), intent(out) :: ierr
  call ira_svd_x( nat2, typ1_in, coords1_in, nat2, typ2_in, coords2_in, &
       kmax_factor, rotation, translation, permutation, hd, rmsd, ierr )
end subroutine ira_svd

subroutine cshda_svd( nat1, typ1_in, coords1_in, &
     nat2, typ2_in, coords2_in, &
     dthr, recenter, &
     perm, dists, rmat, tr, ierr )
  use ira_precision
  use ira_mod, only: cshda_svd_x=>cshda_svd
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), intent(in) :: typ1_in(nat1)
  real(rp),    intent(in) :: coords1_in(3,nat1)
  integer(ip), intent(in) :: nat2
  integer(ip), intent(in) :: typ2_in(nat2)
  real(rp),    intent(in) :: coords2_in(3,nat2)
  real(rp),    intent(in) :: dthr
  logical,     intent(in) :: recenter
  integer(ip), intent(out) :: perm(nat2)
  real(rp),    intent(out) :: dists(nat2)
  real(rp),    intent(out) :: rmat(3,3)
  real(rp),    intent(out) :: tr(3)
  integer(ip), intent(out) :: ierr
  call cshda_svd_x( nat1, typ1_in, coords1_in, nat2, typ2_in, coords2_in, &
       dthr, recenter, perm, dists, rmat, tr, ierr )
end subroutine cshda_svd

subroutine svdrot_m(nat1, typ1, coords1_in, &
     nat2, typ2, coords2_in, &
     rmat, translate, ierr)
  use ira_precision
  use ira_mod, only: svdrot_m_x=>svdrot_m
  implicit none
  integer(ip), intent(in) :: nat1
  integer(ip), dimension(nat1), intent(in) :: typ1
  real(rp), dimension(3, nat1), intent(in) :: coords1_in
  integer(ip), intent(in) :: nat2
  integer(ip), dimension(nat2), intent(in) :: typ2
  real(rp), dimension(3, nat2), intent(in) :: coords2_in
  real(rp), dimension(3, 3), intent(out) :: rmat
  real(rp), dimension(3), intent(out) :: translate
  integer(ip), intent(out) :: ierr
  call svdrot_m_x( nat1, typ1, coords1_in, nat2, typ2, coords2_in, &
       rmat, translate, ierr )
end subroutine svdrot_m

subroutine ira_get_err_msg(ierr, msg)
  use ira_precision
  use ira_mod, only: ira_get_err_msg_x=>ira_get_err_msg
  integer(ip), intent(in) :: ierr
  character(512), intent(out) :: msg

  call ira_get_err_msg_x(ierr, msg)
end subroutine ira_get_err_msg
