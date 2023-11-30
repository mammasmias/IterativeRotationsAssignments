program main
  implicit none

  integer :: nat1, nat2, dnat, nat_pass
  integer, allocatable :: typ1(:), typ2(:), typ_pass(:)
  real, allocatable :: coords1(:,:), coords2(:,:), coords_pass(:,:)
  integer :: i, ierr
  integer, allocatable :: candidate_1(:), candidate_2(:)
  real, dimension(3) :: translation, svd_tr
  real, dimension(3,3) :: rotation, svd_rot
  real :: hd_out, hd
  integer, allocatable :: permutation(:)
  real :: kmax_factor
  real :: rmsd
  real, dimension(3) :: rdum

  !!
  !! remember! structure 1 =< structure 2
  !!


  !!
  !! read xyz structure 1
  !!
  read(*,*) nat1
  allocate( typ1(1:nat1))
  allocate( coords1(1:3,1:nat1))
  read(*,*)
  do i = 1, nat1
     read(*,*) typ1(i),coords1(1,i),coords1(2,i),coords1(3,i)
  end do

  !!
  !! read xyz structure 2
  !!
  read(*,*) nat2
  allocate( typ2(1:nat2))
  allocate( coords2(1:3,1:nat2))
  read(*,*)
  do i = 1, nat2
     read(*,*) typ2(i),coords2(1,i),coords2(2,i),coords2(3,i)
  end do

  !!
  !! make sure that struc1 <= struc2, otherwise switch
  !!
  if( nat1 .gt. nat2 ) then
     !! switch 1 to pass
     nat_pass = nat1
     call move_alloc( typ1, typ_pass )
     call move_alloc( coords1, coords_pass )
     !! switch 2 to 1
     nat1 = nat2
     call move_alloc( typ2, typ1 )
     call move_alloc( coords2, coords1 )
     !! switch pass to 2
     nat2 = nat_pass
     call move_alloc( typ_pass, typ2 )
     call move_alloc( coords_pass, coords2 )
  endif





  !!
  !! =================================
  !! The following routine ira_svd encapsulates a default call to IRA. This includes
  !! setting the candidate central atoms as per default.
  !! To see a deconstructed sequence of calls, see the commented code below.
  !! =================================
  !!
  allocate( permutation(1:nat2))
  kmax_factor = 1.9
  call ira_svd( nat1, typ1, coords1, &
                nat2, typ2, coords2, &
                kmax_factor, rotation, translation, permutation, hd, rmsd, ierr )

  if( ierr /= 0 ) return

  !!
  !! apply to structure 2
  !!
  typ2(:) = typ2(permutation(:))
  coords2(:,:) = coords2(:,permutation(:))
  do i = 1, nat2
     coords2(:,i) = matmul(rotation,coords2(:,i)) + translation
  end do
  !!
  !! or alternatively apply to structure 1
  !!
  ! do i = 1, nat1
  !    coords1(:,i) = matmul(transpose(rotation),coords1(:,i)) - matmul(transpose(rotation),translation)
  ! end do

  write(*,*) nat1
  write(*,*)
  do i = 1, nat1
     write(*,*) typ1(i), coords1(:,i)
  end do
  write(*,*) nat2
  write(*,*)
  do i = 1, nat2
     write(*,*) typ2(i),coords2(:,i)
  end do

  hd = 0.0
  do i = 1, nat1
     hd = max( hd, norm2(coords1(:,i) - coords2(:,i)) )
  end do

  rmsd = 0.0
  do i = 1, nat1
     rdum = coords1(:,i) - coords2(:,i)
     rmsd = rmsd + dot_product(rdum,rdum)
  end do

  write(*,*) 'final values of dH and rmsd:'
  write(*,*) hd, sqrt(rmsd/nat1)

  write(*,*)
  write(*,*) "rotation matrix:"
  write(*,'(3f9.4)') rotation(1,:)
  write(*,'(3f9.4)') rotation(2,:)
  write(*,'(3f9.4)') rotation(3,:)

  write(*,*)
  write(*,*) "translation vector"
  write(*,*) translation



  !!
  !! =================================
  !! The following code shows a deconstructed series of calls
  !! to the main algorithm routines.
  !! Attention, use with care and try to understand what you are doing
  !! before applying them.
  !! =================================
  !!



  ! !!
  ! !! form candidates for central atms in 1 and 2
  ! !!
  ! allocate( candidate_1(1:nat1) )
  ! allocate( candidate_2(1:nat2) )
  ! call set_candidates( nat1, typ1, coords1, &
  !                      nat2, typ2, coords2, &
  !                      candidate_1, candidate_2 )

  ! !!
  ! !! call main ira with candidates
  ! !!
  ! allocate( permutation(1:nat2))
  ! kmax_factor = 1.2
  ! call ira_unify( nat1, typ1, coords1, candidate_1, &
  !      nat2, typ2, coords2, candidate_2, &
  !      kmax_factor, rotation, translation, permutation, hd_out)

  ! !!
  ! !! apply found transformation
  ! !!
  ! typ2(:) = typ2(permutation(:))
  ! coords2(:,:) = coords2(:,permutation(:))
  ! !!
  ! do i = 1, nat2
  !    coords2(:,i) = matmul(rotation, coords2(:,i)) + translation
  ! end do

  ! !!
  ! !! call SVD
  ! !!
  ! call svdrot_m( nat1, typ1, coords1, &
  !      nat1, typ2(1:nat1), coords2(:,1:nat1), &
  !      svd_rot, svd_tr )

  ! !!
  ! !! apply svd
  ! !!
  ! ! do i = 1, nat2
  ! !    coords2(:,i) = matmul( svd_rot, coords2(:,i) ) + svd_tr
  ! ! end do

  ! !! rotate back to orig
  ! do i = 1, nat2
  !    coords2(:,i) = matmul(transpose(rotation),coords2(:,i)) - matmul(transpose(rotation),translation)
  ! end do

  ! !!
  ! !! put together apx and svd
  ! !!
  ! rotation = matmul(svd_rot, rotation)
  ! translation = matmul(svd_rot,translation) + svd_tr
  ! !! apply
  ! do i = 1, nat2
  !    coords2(:,i) = matmul(rotation,coords2(:,i)) + translation
  ! end do
  ! ! do i = 1, nat1
  ! !    coords1(:,i) = matmul(transpose(rotation),coords1(:,i)) - matmul(transpose(rotation),translation)
  ! ! end do

end program main
