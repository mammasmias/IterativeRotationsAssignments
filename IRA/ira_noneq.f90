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

program ira_noneq
  !!
  !! IRA algorithm for an input which contains different number of points.
  !! The central atom for structure 1 is given in c_atm1.
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

  !!
  if( nat2 .lt. nat1 ) then
     write(*,*) 'WARNING:'
     write(*,*) '      Structure 2 has less atoms than structure 1.'
     write(*,*) '      Try swapping them. Stopping.'
     stop
  endif


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
  !! different number of atoms, call ira with central atom in structure 1,
  !! with index c_atm1.
  !! When a better guess is possible, this can be modified.
  !!
  c_atm1 = 1
  !!
  call ira_nonequal( nat1, typ1, coords1, &
                     nat2, typ2, coords2, &
                     c_atm1, apx_rot, apx_tr, apx_p )

  !!
  !! apply the found permutation
  !!
  typ2(:) = typ2( apx_p(:) )
  coords2(:,:) = coords2(:, apx_p(:) )
  !!
  !! apply the found apx transform (not strictly necessary)
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

end program ira_noneq


