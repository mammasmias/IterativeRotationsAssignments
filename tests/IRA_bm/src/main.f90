#include "dbg.h"
program main
  use ira_mod, only: ira_rp, ira_svd, sofi_construct_operation
  use dbg
  use routines
  implicit none
  integer :: nat1, nat2
  character(len=4), allocatable :: ctyp1(:), ctyp2(:)
  integer, allocatable :: ityp1(:), ityp2(:)
  real(ira_rp), allocatable :: pos1(:,:), pos2(:,:)

  type( t_settings ) :: settings
  integer :: i, ifile, ierr, itry, ntry
  real(ira_rp) :: dt0, dt1, dtfile
  real(ira_rp) :: rmat(3,3), tr(3), expected_rmat(3,3), expected_tr(3), &
       rmat_out(3,3), tr_out(3)
  integer, allocatable :: perm(:), expected_perm(:), perm_out(:)
  real(ira_rp) :: hd, rmsd, time_tot
  real(ira_rp) :: ax(3), an, z, gc(3)
  type( t_sym ) :: sym

  ! same random sequence on every run, so failures are reproducible
  call random_init( repeatable=.true., image_distinct=.true. )

  ntry = 50
  time_tot = 0.0_ira_rp
  call settings%parse_args()
  do ifile = 1, settings%nfiles
     !
     dtfile = 0.0_ira_rp
     ! read struc from file
     ierr = read_xyz( trim(settings%fnames(ifile)), nat1, ctyp1, pos1 )
     if( ierr /= 0 ) then
        write(*,*) "problem with file: "//trim(settings%fnames(ifile))
        error stop 1
     end if
     ! copy to struc2
     nat2 = nat1
     allocate( ityp1(nat1) )
     allocate( ityp2(nat2) )
     allocate( ctyp2, source=ctyp1 )
     allocate( pos2, source=pos1 )
     allocate( perm(1:nat1))
     allocate( perm_out(1:nat1) )
     ! convert ctyp to ityp
     call str_to_int( ctyp1, ctyp2, ityp1, ityp2 )

     ! compute sym, to get degenerate solutions
     sym = sofi_compute( nat1, ityp1, pos1, 0.2_ira_rp, ierr )
     ASSERT( ierr == 0, "sofi_compute failed, ierr="//tostr(ierr) )
     ! sofi works on the gc-centered structure
     gc = sum( pos1, 2 )/real( nat1, ira_rp )

     do itry = 1, ntry

        ! generate random transformations:
        ! > permutation
        call random_permutation( nat1, perm )
        pos2 = pos1(:,perm)
        ityp2 = ityp1(perm)

        ! > rotation/inversion
        call random_number( ax )
        ax = 2.0_ira_rp*ax - 1.0_ira_rp   ! axis in any direction, not only the +++ octant
        call random_number( an )
        call sofi_construct_operation( "C", ax, an, rmat, ierr )
        call random_number(z)
        if( z > 0.5_ira_rp ) rmat(:,3)=-rmat(:,3)
        ! > translation
        call random_number(tr)
        tr = tr*10.0

        ! transform
        do i = 1, nat2
           pos2(:,i) = matmul( rmat, pos2(:,i) ) + tr
        end do



        ! save info
        expected_rmat = transpose(rmat)
        expected_tr = -matmul(transpose(rmat), tr)
        expected_perm = inverse_perm( perm )

        call cpu_time(dt0)
        ! call match
        call ira_svd( nat1, ityp1, pos1, &
             nat2, ityp2, pos2, &
             1.2_ira_rp, rmat_out, tr_out, perm_out, hd, rmsd, ierr )
        call cpu_time(dt1)

        ! check outcome
        if( settings%log_each ) then
           write(*,*) "> itry:",itry, "ierr=",ierr
           write(*,*) "nbas:", sym%nbas
           write(*,*) "hd, rmsd:",hd, rmsd
           write(*,*) "rotation matrix;         expected:"
           write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") rmat_out(1,:), expected_rmat(1,:)
           write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") rmat_out(2,:), expected_rmat(2,:)
           write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") rmat_out(3,:), expected_rmat(3,:)
           write(*,*) "translation vector; expected:"
           write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") tr_out, expected_tr
           write(*,*) "permutation"
           write(*,"(10(i3,1x))") perm_out
           write(*,*) "permutation expected"
           write(*,"(10(i3,1x))") expected_perm
        end if

        ! assert
        ASSERT( ierr == 0, "ira_svd failed, ierr="//tostr(ierr) )
        !! expecting perfect match always
        ASSERT( within(hd, 0.0_ira_rp, 1e-5), "unexpected hd:"//tostr(hd))
        ASSERT( within(rmsd, 0.0_ira_rp, 1e-5), "unexpected rmsd:"//tostr(rmsd))
        !! testing on rmat, tr, and perm can be ambiguous if struc has symmetry
        ASSERT( check_permutation(perm_out)==0, "permutation is not a bijection:"//tostr(perm_out,"(*(i3,1x))") )

        block
          ! With symmetric strucs, any solution R' = S * R^T is valid, where S is a
          ! symmetry op of struc1 (S=E gives the plain expected values).
          ! Sofi convention:  S c(:,perm_list(j)) = c(:,j),  with c = pos1 - gc
          ! which gives:
          !   rmat = S * expected_rmat
          !   perm = expected_perm( perm_list(:,isym) )
          !   tr   = S * (expected_tr - gc) + gc
          integer :: isym
          logical :: found
          real(ira_rp) :: s(3,3), rmat_sym(3,3), tr_sym(3)
          integer, allocatable :: perm_sym(:)
          found = .false.
          do isym = 1, sym%nbas
             s = sym%bas_list(:,:,isym)
             rmat_sym = matmul( s, expected_rmat )
             if( .not. all(within(rmat_out, rmat_sym, 1e-5)) ) cycle
             found = .true.
             perm_sym = expected_perm( sym%perm_list(1:nat1,isym) )
             tr_sym   = matmul( s, expected_tr - gc ) + gc
             ASSERT( all(within(perm_out, perm_sym, 0.0)), "isym="//tostr(isym)//" expected perm:"//tostr(perm_sym,"(*(i3,1x))")//new_line("a")//"> got          :"//tostr(perm_out,"(*(i3,1x))"))
             ASSERT( all(within(tr_out, tr_sym, 1e-5)), "isym="//tostr(isym)//" expected tr:"//tostr(tr_sym, "(3(g0.6,1x))")//new_line("a")//"> got        :"//tostr(tr_out,"(3(g0.6,1x))") )
             exit
          end do
          ASSERT( found, "rmat_out is not expected_rmat times any symmetry op, nbas="//tostr(sym%nbas)//new_line("a")//"> got: "//tostr(reshape(rmat_out,[9]),fmt="(9(g0.6,1x))") )
        end block


        dtfile = dtfile + (dt1-dt0)

     end do
     deallocate( ctyp1, ityp1, pos1 )
     deallocate( ctyp2, ityp2, pos2 )
     deallocate( perm, expected_perm )
     deallocate( perm_out )
     time_tot = time_tot + dtfile
     if( settings%log_file ) then
        write(*,"(a,2x, f12.8, 2x, f12.8)") &
             trim(settings%fnames(ifile)), dtfile, dtfile/real(ntry,ira_rp)
     end if
  end do

  write(*,"(a,2x,g0.8)") "time_tot", time_tot

end program main
