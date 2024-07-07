program sofi
  !> @details
  !!
  !! This is a program which returns the symmetry operations
  !! of a structure given in input. It also prints the Point Group associated
  !! with the list of found SymmOps.
  !!
  !!
  implicit none

  integer :: nat                             !< @brief number of atoms
  integer, allocatable :: typ(:)             !< @brief integer of atomic types
  character(len=4), allocatable :: styp(:)   !< @string of atomic types
  real, allocatable :: coords(:,:)           !< @brief atomic positions

  integer :: i, j

  real, dimension(3) :: gc                   !< @brief origin point (geometric center)
  real :: sym_thr                            !< @brief threshold for symmetry operations
  real, allocatable :: bas_list(:,:,:)       !< @brief list of symmetry operations
  integer, allocatable :: perm_list(:,:)     !< @brief list of permutations for each symmetry operation
  integer :: nbas                            !< @brief total number of symmetry operations
  integer :: nmax                            !< @brief max space for allocation (see sofi_tools.f90)
  character(len=10) :: pg                    !< @brief the point group
  integer :: n_prin_ax                       !< @brief number of principal axes
  real, allocatable :: prin_ax(:,:)          !< @brief list of principal axes
  logical :: prescreen_ih

  real, allocatable :: angle_out(:)          !< @brief list of angles of symm operations
  real, allocatable :: ax_out(:,:)           !< @brief list of axes of symmetry operations
  character(len=1), allocatable :: op_out(:) !< @brief list of Op
  integer, allocatable :: n_out(:)           !< @brief list of n values
  integer, allocatable :: p_out(:)           !< @brief list of p values
  real, allocatable :: dHausdorff_out(:)           !< @brief list of dHausdorff values
  integer :: ierr
  character(128) :: msg

  !! nmax is imposed from sofi_tools.f90
  nmax = 400

  !! threshold for sym
  sym_thr = 0.05

  read(*,*) nat
  read(*,*)
  allocate( typ(1:nat), source=0)
  allocate( coords(1:3,1:nat), source=0.0)
  allocate( styp(1:nat), source="aaaa")
  do i = 1, nat
     !! read at typ into string
     read(*,*) styp(i), coords(:,i)
  end do
  !! convert atomic typ from string into integer
  call str_to_int( styp, typ )

  !! recenter to gc
  gc = sum( coords(:,:),2)/nat
  do i = 1, nat
     coords(:,i) = coords(:,i) - gc
  end do


  write(*,*) nat
  write(*,*)
  do i = 1, nat
     write(*,*) styp(i), coords(:,i), typ(i)
  end do

  !! allocate to nmax just for space. The routine get_symmops expects size 1:nmax,
  !! all entries except for the first nbas elements are zero on output
  allocate( bas_list(1:3, 1:3, 1:nmax))
  allocate( perm_list(1:nat, 1:nmax))
  allocate( op_out(1:nmax))
  allocate( n_out(1:nmax))
  allocate( p_out(1:nmax))
  allocate( ax_out(1:3, 1:nmax))
  allocate( angle_out(1:nmax))
  allocate( dHausdorff_out(1:nmax))
  allocate( prin_ax(1:3, 1:nmax))
  prescreen_ih = .False.

  call sofi_compute_all( nat, typ, coords, sym_thr, prescreen_ih,  &
       nbas, bas_list, perm_list, op_out, n_out, p_out, ax_out, angle_out, dHausdorff_out, pg,&
       n_prin_ax, prin_ax, ierr )
  if( ierr /= 0 ) then
     write(*,*) "f90 prog got nonzero ierr:", ierr
     call sofi_get_err_msg( ierr, msg )
     write(*,*) trim(msg)
     return
  end if


  do i = 1, nbas
     write(*,'(2x,i0)') i
     write(*,'(2x,a,3x,a2,x,i0,"^",i0)') "operation:",op_out(i), n_out(i), p_out(i)
     write(*,'(2x, a,f9.4)') "angle",angle_out(i)
     write(*,'(2x, a,3f9.5)') "axis",ax_out(:,i)
     write(*,'(2x,a)') "matrix:"
     do j = 1, 3
        write(*,'(3f12.6)') bas_list(j,:,i)
     end do
     write(*,'(2x,a,f12.7)') "dHausdorff",dHausdorff_out(i)
     write(*,'(2x, a)') "permutation of atoms:"
     write(*,'(20i4)') perm_list(:,i)
     write(*,*)
  end do

  write(*,*) "PG:",trim(pg)
  write(*,"(a,1x,i0)") "List of principal axes, N_prin_ax =",n_prin_ax
  do i = 1, n_prin_ax
     write(*,"(3f9.4)") prin_ax(:,i)
  end do
  deallocate( bas_list, perm_list, op_out, n_out, p_out, ax_out, angle_out, dHausdorff_out, prin_ax )

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
  deallocate( styp )



contains

  subroutine str_to_int( str_arr, int_arr )

    implicit none
    character(len=*), dimension(:), intent(in) :: str_arr
    integer, dimension(:), intent(out) :: int_arr

    integer :: i, j, n, wlen
    character(:), allocatable :: old(:)
    logical :: isnew
    integer :: nnew

    if( size(int_arr) /= size(str_arr) ) then
       write(*,*) "unequal size arrays in str_to_int. Stopping"
       stop
    end if

    wlen = len(str_arr(1))
    n = size(int_arr, 1)
    allocate( character(len=wlen)::old(n) )

    old(:) = "xcba"
    nnew = 1
    old(nnew) = str_arr(1)
    do i = 2, n
       !! check if element is new
       isnew = .not. any( old .eq. str_arr(i) )
       if( isnew ) then
          nnew = nnew + 1
          old( nnew ) = str_arr(i)
       end if
    end do

    do i = 1, n
       do j = 1, nnew
          if( str_arr(i) .eq. old(j) ) int_arr(i) = j
       end do
    end do


    deallocate(old)
  end subroutine str_to_int


end program sofi

