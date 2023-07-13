program sofi
  !!
  !! This is a program which returns the symmetry operations
  !! of a structure given in input. It also prints the Point Group associated
  !! with the list of found SymmOps.
  !!
  !!
  implicit none

  integer :: nat
  integer, allocatable :: typ(:)
  real, allocatable :: coords(:,:)
  integer :: i, j
  real, dimension(3) :: gc
  real :: sym_thr
  real, allocatable :: bas_list(:,:,:)
  integer, allocatable :: perm_list(:,:)
  integer :: nbas
  integer :: nmax
  character(len=10) :: pg

  real, allocatable :: angle_out(:), ax_out(:,:)
  character(len=2), allocatable :: op_out(:)

  integer, allocatable :: n_out(:), p_out(:)
  real, allocatable :: dmax_out(:)

  !! nmax is imposed from sofi_tools.f90
  nmax = 200

  !! threshold for sym
  sym_thr = 0.02

  read(*,*) nat
  read(*,*)
  allocate( typ(1:nat), source=0)
  allocate( coords(1:3,1:nat), source=0.0)
  do i = 1, nat
     read(*,*) typ(i), coords(:,i)
  end do

  !! recenter to gc
  gc = sum( coords(:,:),2)/nat
  do i = 1, nat
     coords(:,i) = coords(:,i) - gc
  end do


  write(*,*) nat
  write(*,*)
  do i = 1, nat
     write(*,*) typ(i), coords(:,i)
  end do

  !! allocate to nmax just for space. The routine get_symmops expects size 1:nmax,
  !! all entries except for the first nbas are zero on output
  allocate( bas_list(1:3,1:3,1:nmax))
  allocate( perm_list(1:nat, 1:nmax))
  allocate( op_out(1:nmax))
  allocate( n_out(1:nmax))
  allocate( p_out(1:nmax))
  allocate( ax_out(1:3,1:nmax))
  allocate( angle_out(1:nmax))
  allocate( dmax_out(1:nmax))
  call sofi_compute_all( nat, typ, coords, sym_thr, &
       nbas, bas_list, perm_list, op_out, n_out, p_out, ax_out, angle_out, dmax_out, pg )

  do i = 1, nbas
     write(*,'(2x,i0)') i
     write(*,'(2x,a,3x,a2,x,i0,"^",i0)') "operation:",op_out(i), n_out(i), p_out(i)
     write(*,'(2x, a,f9.4)') "angle",angle_out(i)
     write(*,'(2x, a,3f9.4)') "axis",ax_out(:,i)
     write(*,'(2x,a)') "matrix:"
     do j = 1, 3
        write(*,'(3f12.6)') bas_list(j,:,i)
     end do
     write(*,'(2x,a,f9.4)') "dmax",dmax_out(i)
     write(*,'(2x, a)') "permutation of atoms:"
     write(*,'(25i3)') perm_list(:,i)
     write(*,*)

  end do

  write(*,*) "PG:",trim(pg)
  deallocate( bas_list, perm_list, op_out, n_out, p_out, ax_out, angle_out, dmax_out )

  ! allocate( bas_list(1:3,1:3,1:nmax))
  ! allocate( perm_list(1:nat, 1:nmax))

  ! !! find the symmetry operations and associated permutations
  ! call sofi_get_symmops( nat, typ, coords, sym_thr, nbas, bas_list, perm_list )


  ! !! output symmetry operations
  ! write(*,*) "Number of SymmOps found:",nbas
  ! do i = 1, nbas
  !    rmat = bas_list(:,:,i)
  !    call sofi_analmat( rmat, op, n, p, ax, angle )
  !    write(*,*) repeat('=',30)
  !    write(*,'(a9,a4,g0,a1,g0,a8,3f9.4,a8,x,f8.4)') "SymmOp:",op, n, '^', p, "axis:",ax, "angle:", angle
  !    write(*,*) 'associated matrix:'
  !    do j = 1, 3
  !       write(*,'(3f12.6)') rmat(j,:)
  !    end do
  ! end do

  ! write(*,*)
  ! write(*,*)


  ! write(*,*) "number of SymmOps entering get_pg:",nbas

  ! verb = .false.
  ! call sofi_get_pg( nbas, bas_list, pg, verb )
  ! write(*,*) repeat('=',20)
  ! write(*,*) 'PG is: ',pg


  ! !! write the list of unique axes and angles
  ! allocate( op_out(1:nbas))
  ! allocate( angle_out(1:nbas))
  ! allocate( ax_out(3,nbas))
  ! call sofi_unique_ax_angle( nbas, bas_list, op_out, ax_out, angle_out )
  ! do i = 1, nbas
  !    write(*,'(i4,x,a,x,3f9.4,4x,f9.4)') i, op_out(i), ax_out(:,i), angle_out(i)
  ! end do
  ! deallocate( op_out, angle_out, ax_out )


  ! deallocate( bas_list, perm_list )



  deallocate( typ, coords )

end program sofi
