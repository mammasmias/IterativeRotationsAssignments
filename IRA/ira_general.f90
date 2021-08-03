!!
!! Copyright (C) 2021, MAMMASMIAS Consortium
!! Written by: Miha Gunde
!!
!! This file is distributed under the terms of GNU General Public License.
!! See the license at: http://www.gnu.org/licenses/gpl-3.0.txt
!!

program ira_general
  !!
  !! IRA algorithm for an input which can contain equal or different number of points.
  !! When equal number of atoms, the origin is given by the geometrical center, and
  !! when different number of atoms, the central atom for structure 1 is given in c_atm1.
  !!
  !! =========================
  !! nat1     -> number of atoms in structure 1
  !! nat2     -> number of atoms in structure 2
  !! typ1     -> atomic types of structure 1, stored as integers
  !! typ2     -> atomic types of structure 2, stored as integers
  !! coords1  -> atomic coordinates of structure 1
  !! coords2  -> atomic coordinates of structure 2
  !! apx_rot  -> 3x3 rotation matrix from IRA
  !! apx_tr   -> 3dim translation vector from IRA
  !! svd_rot  -> 3x3 rotation matrix from SVD
  !! svd_tr   -> 3dim translation vector from SVD
  !! apx_p    -> array of permutation from IRA
  use read_typ
  implicit none
  integer :: nat1, nat2
  integer, allocatable :: typ1(:), typ2(:)
  real, allocatable :: coords1(:,:), coords2(:,:)
  integer :: i, c_atm1
  real, dimension(3,3) :: apx_rot, svd_rot
  real, dimension(3) :: apx_tr, svd_tr, rdum
  integer, allocatable :: apx_p(:)
  real :: rmsd
  character(len=2) :: ctyp

  !!
  !! read xyz structure 1
  !!
  read(*,*) nat1
  allocate( typ1(1:nat1))
  allocate( coords1(1:3,1:nat1))
  read(*,*)
  do i = 1, nat1
     read(*,*) ctyp,coords1(1,i),coords1(2,i),coords1(3,i)
     typ1(i) = get_int_typ( ctyp )
  end do

  !!
  !! read xyz structure 2
  !!
  read(*,*) nat2
  allocate( typ2(1:nat2))
  allocate( coords2(1:3,1:nat2))
  read(*,*)
  do i = 1, nat2
     read(*,*) ctyp,coords2(1,i),coords2(2,i),coords2(3,i)
     typ2(i) = get_int_typ( ctyp )
  end do


  !!==========================================================
  !! The IRA procedure starts here
  !!
  !!
  !! initial P = diagonal
  allocate( apx_p(1:nat2) )
  do i = 1, nat2
     apx_p(i) = i
  end do
  !!
  !! decide what to do
  !!
  if( nat1 .eq. nat2 ) then
     !!
     !! equal number of atoms, call ira with geo center
     !!
     call ira_equal( nat1, typ1, coords1, &
          nat2, typ2, coords2, &
          apx_rot, apx_tr, apx_p )
     !!
  elseif( nat1 .lt. nat2 ) then
     !!
     !! different number of atoms, call ira with central atom in structure 1,
     !! with index c_atm1 (this can be modified).
     !!
     c_atm1 = 1
     !!
     call ira_nonequal( nat1, typ1, coords1, &
                        nat2, typ2, coords2, &
                        c_atm1, apx_rot, apx_tr, apx_p )
     !!
  else
     !!
     write(*,*) 'WARNING:'
     write(*,*) '    Number of atoms in structure 1 cannot be larger than'
     write(*,*) '    number of atoms in structure 2. You can swap the'
     write(*,*) '    structures. Stopping.'
     stop
     !!
  endif

  !!
  !! apply the found permutation
  !!
  typ2(:) = typ2( apx_p(:) )
  coords2(:,:) = coords2(:, apx_p(:) )
  !!
  !! apply the found apx transform (not necessary)
  !!
  ! do i = 1, nat2
  !    coords2(:,i) = matmul(apx_rot,coords2(:,i)) + apx_tr
  ! end do

  !!
  !! get svd for points up to nat1
  !!
  call svdrot_m( nat1, typ1, coords1, &
                 nat1, typ2(1:nat1), coords2(1:3,1:nat1), &
                 svd_rot, svd_tr )

  ! write(*,*) 'svd rotation matrix:'
  ! do i = 1, 3
  !    write(*,*) svd_rot(i,:)
  ! end do

  !!
  !! apply the transformation found by svd
  !!
  do i = 1, nat2
     coords2(:,i) = matmul( svd_rot, coords2(:,i) ) + svd_tr
  end do
  !!
  !!
  !! This is the end of IRA procedure
  !!==========================================================



  !!
  !! output matched structures
  !!
  write(*,*) nat1
  write(*,*) 'structure 1 as input'
  do i = 1, nat1
     write(*,*) all_ctyp( typ1(i) ),coords1(:,i)
  end do

  write(*,*) nat2
  write(*,*) 'structure 2 matched'
  do i = 1, nat2
     write(*,*) all_ctyp( typ2(i) ),coords2(:,i)
  end do

  !! get final rmsd up to nat1
  rmsd = 0.0
  do i = 1, nat1
     rdum = coords1(:,i) - coords2(:,i)
     rmsd = rmsd + dot_product(rdum,rdum)
  end do
  write(*,'(a,f6.3)') 'final rmsd:',sqrt(rmsd/nat1)

  deallocate( typ1, coords1 )
  deallocate( typ2, coords2 )
  deallocate( apx_p )

end program ira_general


