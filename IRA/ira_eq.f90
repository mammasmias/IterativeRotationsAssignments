!!
!! Copyright (C) 2021, MAMMASMIAS Consortium
!! Written by: Miha Gunde
!!
!! SPDX-License-Identifier: GPL-3.0-or-later
!! SPDX-License-Identifier: Apache-2.0
!! See the file LICENSE.txt for further information.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!

program ira_eq
  !!
  !! IRA match two structures containing the same number of points.
  !! The center is set on the geometrical center of both structures.
  !!
  use read_typ
  implicit none
  integer :: nat1, nat2
  integer, allocatable :: typ1(:), typ2(:)
  real, allocatable :: coords1(:,:), coords2(:,:)
  integer :: i
  real, dimension(3,3) :: apx_rot, svd_rot
  real, dimension(3) :: apx_tr, svd_tr, rdum
  integer, allocatable :: apx_p(:)
  real :: rmsd
  character(len=2) :: ctyp

  !! read xyz structure 1
  read(*,*) nat1
  allocate( typ1(1:nat1))
  allocate( coords1(1:3,1:nat1))
  read(*,*)
  do i = 1, nat1
     read(*,*) ctyp,coords1(1,i),coords1(2,i),coords1(3,i)
     typ1(i) = get_int_typ( ctyp )
  end do

  !! read xyz structure 2
  read(*,*) nat2
  allocate( typ2(1:nat2))
  allocate( coords2(1:3,1:nat2))
  read(*,*)
  do i = 1, nat2
     read(*,*) ctyp,coords2(1,i),coords2(2,i),coords2(3,i)
     typ2(i) = get_int_typ( ctyp )
  end do

  !!
  if( nat1 .ne. nat2 ) then
     write(*,*) 'WARNING:'
     write(*,*) '     Nonequal number of atoms in structures:',nat1, nat2
     write(*,*) '     Stopping.'
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
  !! ira with equal number of atoms, origin at geoemtric center
  !!
  call ira_equal( nat1, typ1, coords1, &
                  nat2, typ2, coords2, &
                  apx_rot, apx_tr, apx_p )

  !!
  !! apply found permutation
  !!
  ! typ2(:) = typ2( apx_p(:) )
  ! coords2(:,:) = coords2(:, apx_p(:) )
  call permute_int_1d( nat2, typ2, apx_p )
  call permute_real_2d( nat2, 3, coords2, apx_p )
  !!
  !! apply apx transformation (not strictly necessary)
  !!
  do i = 1, nat2
     coords2(:,i) = matmul(apx_rot,coords2(:,i)) + apx_tr
  end do

  !!
  !! get svd
  !!
  call svdrot_m( nat1, typ1, coords1, &
                 nat2, typ2, coords2, &
                 svd_rot, svd_tr )

  ! write(*,*) 'svd rotation matrix:'
  ! do i = 1, 3
  !    write(*,*) svd_rot(i,:)
  ! end do

  !!
  !! apply transformation found by SVD
  !!
  do i = 1, nat2
     coords2(:,i) = matmul( svd_rot, coords2(:,i) ) + svd_tr
  end do
  !!
  !!
  !! This is the end of IRA procedure
  !!==========================================================


  !! output matched structures
  write(*,*) nat1
  write(*,*) 'conf1 as input'
  do i = 1, nat1
     write(*,*) typ1(i),coords1(:,i)
  end do

  write(*,*) nat2
  write(*,*) 'conf2 matched'
  do i = 1, nat2
     write(*,*) typ2(i),coords2(:,i)
  end do

  !! get final rmsd
  rmsd = 0.0
  do i = 1, nat1
     rdum = coords1(:,i) - coords2(:,i)
     rmsd = rmsd + dot_product(rdum,rdum)
  end do
  write(*,'(a,f6.3)') 'final rmsd:',sqrt(rmsd/nat1)

  deallocate( typ1, coords1 )
  deallocate( typ2, coords2 )
  deallocate( apx_p )

end program ira_eq


