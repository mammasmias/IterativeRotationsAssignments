program pgsimple

  integer :: nat
  integer, allocatable :: typ(:)
  real, allocatable :: coords(:,:)

  real :: sym_thr
  character(len=10) :: pg
  integer :: i
  logical :: verb

  read(*,*) nat

  allocate( typ(1:nat))
  allocate( coords(1:3,1:nat))

  !! read the symmetry threshold
  ! read(*,*) sym_thr
  read(*,*)

  !! hard-code the symmetry threshold
  sym_thr = 0.05

  do i = 1, nat
     read(*,*) typ(i), coords(:,i)
  end do

  !! control verbosity of the routine
  verb = .false.
  ! verb = .true.

  !! call the sofi_struc_pg routine, which does the shift to geometric center
  !! inside the routine, no need to do it here.
  call sofi_struc_pg( nat, typ, coords, sym_thr, pg, verb )

  write(*,*) 'PG is:', pg

  deallocate( typ, coords )
end program pgsimple
