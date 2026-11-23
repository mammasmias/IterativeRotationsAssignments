#include "dbg.h"
program test_ira
  use dbg
  ! use ira_precision, only: ira_rp => rp
  use ira_mod

  implicit none

  integer :: nat1, nat2


  nat1 = 13; nat2 = 13
  write(*,*) repeat("=",60)
  write(*,*) "Test match on randomized structures with nat="//tostr(nat1)
  write(*,*) "Expecting perfect match."
  write(*,*) repeat("=",60)
  block
    integer :: typ1(nat1), typ2(nat2)
    real(ira_rp) :: coords1(3,nat1), coords2(3,nat2)
    real(ira_rp) :: rmat(3,3), rmat_expected(3,3)
    real(ira_rp) :: tr(3), tr_expected(3)
    real(ira_rp) :: hd, hd_expected, rmsd, rmsd_expected
    integer :: perm(nat2), perm_expected(nat2)
    integer :: ierr
    integer :: itry, ntry, i
    real(ira_rp) :: an, ax(3), z
    typ1 = 1; typ2 = typ1
    ntry = 50

    coords1 = real( reshape([ &
         -0.49580341431341102     ,  0.97081809720047896     ,  0.37341428393206699,     &
         -1.0561165586604031     , -0.47245030034569702     , -0.37449783746929599,     &
        -0.63509643604830102     , -0.66670775781629399     ,  0.66219897224722502,     &
        -0.83642178463603800     ,  0.59155936354314198     , -0.64507703368827696,     &
         0.59636158974191400     ,  0.80558700835169195     ,  0.23843961936244301,     &
         0.25975283602565602     ,  0.71540296586945495     , -0.78971023572037002,     &
         -9.7433080485659004E-002,  -1.0381280426042809     , -0.31233048725326901,     &
          9.2545024674446003E-002,  0.20016737720269201     ,   1.0302106792521439,     &
        -0.18424967464368000     , -0.24756756715313999     ,  -1.0721752232720581,     &
         0.46705990530777303     , -0.73516435209754005     ,  0.56288324725063099,     &
         0.82832153270028297     , -0.18917103116725101     , -0.36543024542442698,     &
          1.2108515516656180     ,   4.1719387201624998E-002,  0.70259140999599501,     &
        -0.14977149132819800     ,   2.3934851915116000E-002,  -1.0517149112811999E-002], &
        [3,nat1]), ira_rp )

    do itry = 1, ntry
       ! randomly permute coords1
       call random_permutation(nat1, perm)
       coords2 = coords1(:,perm)
       typ2 = typ1(perm)
       perm_expected = inverse_perm( perm )

       ! randomly rotate
       call random_number( ax )
       ax = ax*2.0 - [1.0, 1.0, 1.0]
       call random_number(an)
       call sofi_construct_operation( "C", ax, an, rmat, ierr )
       ! generate rotoreflection randomly
       call random_number(z)
       if( z > 0.5_ira_rp ) rmat(3,:) = -rmat(3,:)
       rmat_expected = transpose(rmat)


       call random_number(tr)
       tr = tr*10.0
       tr_expected = -matmul(transpose(rmat), tr)

       !! transform
       do i = 1, nat2
          coords2(:,i) = matmul( rmat, coords2(:,i) ) + tr
       end do

       ! compute match
       call ira_svd( nat1, typ1, coords1, &
                     nat2, typ2, coords2, &
                     1.2_ira_rp, rmat, tr, perm, hd, rmsd, ierr )

       write(*,*) "> itry:",itry, "ierr=",ierr
       write(*,*) "hd, rmsd:",hd, rmsd
       write(*,*) "rotation matrix;         expected:"
       write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") rmat(1,:), rmat_expected(1,:)
       write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") rmat(2,:), rmat_expected(2,:)
       write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") rmat(3,:), rmat_expected(3,:)
       write(*,*) "translation vector; expected:"
       write(*,"(3(f12.6,1x),3x,'|',3x,3(f12.6,1x))") tr, tr_expected
       write(*,*) "permutation"
       write(*,"(10(i3,1x))") perm

       !! expecting perfect match always
       ASSERT( within(hd, 0.0_ira_rp, 1e-5), "unexpected hd:"//tostr(hd))
       ASSERT( within(rmsd, 0.0_ira_rp, 1e-5), "unexpected rmsd:"//tostr(rmsd))
       ASSERT( all(within(rmat, rmat_expected, 1e-5)), "expected rmat:"//tostr(reshape(rmat,[9]),fmt="(9(g0.6,1x))")//new_line("a")//"> got          :"//tostr(reshape(rmat_expected,[9]),fmt="(9(g0.6,1x))") )
       ASSERT( all(within(tr, tr_expected, 1e-5)), "expected tr:"//tostr(tr_expected, "(3(g0.6,1x))")//new_line("a")//"> got        :"//tostr(tr,"(3(g0.6,1x))") )
       ASSERT( all(within(perm, perm_expected, 0.0)), "expected perm:"//tostr(perm_expected,"(*(i3,1x))")//new_line("a")//"> got          :"//tostr(perm,"(*(i3,1x))"))
       ASSERT( check_permutation(perm)==0, "permutation is not a bijection:"//tostr(perm,"(*(i3,1x))") )

    end do

  end block

end program test_ira
