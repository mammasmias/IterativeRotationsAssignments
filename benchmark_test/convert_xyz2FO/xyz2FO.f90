program xyz2FO
  !!
  !! small code to transform .xyz format into the data format needed for the FastOverlap
  !! script called "sphericalAlignment.py"
  !!
  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  character(len=2) :: ctyp
  real, allocatable :: coords(:,:)
  integer :: i

  !! read xyz format
  read(*,*) nat
  allocate( typ(1:nat) )
  allocate( coords(1:3, 1:nat) )
  read(*,*)
  do i = 1, nat
     !read(*,*) typ(i), coords(1,i), coords(2,i), coords(3,i)
     read(*,*) ctyp, coords(1,i), coords(2,i), coords(3,i)
  end do

  !! output fastoverlap format
  do i = 1, nat
     write(*,'(3(g0.8,x))') coords(:,i)
  end do

  deallocate( typ, coords )

end program xyz2FO

