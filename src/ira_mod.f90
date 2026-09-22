module ira_mod

  ! the catch-all/facade module

  use ira_precision,  only: &
       ira_ip=>ip,          &
       ira_rp=>rp

  use m_ira_version,  only: &
       ira_get_version

  use m_cshda,        only: &
       cshda,               &
       cshda_from_cost,     &
       cshda_pbc

  use m_ira_routines, only: &
       ira_unify,           &
       ira_svd,             &
       cshda_svd,           &
       svdrot_m,            &
       ira_get_err_msg

  use m_sofi_tools,    only: &
       sofi_nmax => nmax

  use m_sofi_routines, only: &
       sofi_compute_all,     &
       sofi_struc_pg,        &
       sofi_get_symmops,     &
       sofi_get_perm,        &
       sofi_get_combos,      &
       sofi_get_pg,          &
       sofi_analmat,         &
       sofi_ext_Bfield,      &
       sofi_construct_operation, &
       sofi_mat_combos,      &
       sofi_check_collinear, &
       sofi_get_err_msg

  ! private
  public

contains

end module ira_mod
