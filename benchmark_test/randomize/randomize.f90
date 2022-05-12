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
!!===============================================================
!! Description:
!!  This program reads an input structure from a .xyz input, generates a 
!!  random rigid transformation, and returns the input structure transformed
!!  by the generated transformation.
program randomiz
  use read_typ
  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real, allocatable :: coords(:,:)
  integer :: i
  real, dimension(3) :: gc
  character( len = 2) :: ctyp

  !! read xyz input
  read(*,*) nat
  allocate( coords(1:3,1:nat))
  allocate( typ(1:nat))
  read(*,*)
  do i = 1, nat
     read(*,*) ctyp, coords(1,i), coords(2,i),coords(3,i)
     typ(i) = get_int_typ( ctyp )
  end do

  !! take out geo center
  gc = sum(coords(:,:),2)/nat
  do i = 1, nat
      coords(:,i) = coords(:,i) - gc
  end do

  !! randomize
  call randomize( nat, typ, coords )

  !! write randomized structure
  write(*,*) nat
  write(*,*)
  do i = 1, nat
     write(*,*) all_ctyp(typ(i)), coords(:,i)
!     write(*,*) coords(:,i)
  end do


end program randomiz


  subroutine randomize( nat, typ, coords )
    !! generate random rigid transformation: 
    !!      permutation, rotation, translation, reflection,
    !! and apply it to input structure.
    implicit none
    integer, intent(in) :: nat
    integer,dimension(nat), intent(inout) :: typ
    real, dimension(3,nat), intent(inout) :: coords

    integer :: i
    real :: z
    real, dimension(3) :: ax, tr
    integer, dimension(nat) :: p
    real, dimension(3,3) :: rmat

    call set_random_seed()

    do i = 1, nat
       p(i) = i
    end do

    !!
    !! random permute
    !!
    call random_permutation(nat,p)
    ! write(877,*) 'random permutation'
    ! write(877,'(20i4)') p(:)
    typ(:) = typ(p)
    coords(:,:) = coords(:,p)

    !!
    !! random rotate
    !!
    !! random axis
    do i = 1, 3
       call random_number(z)
       z = z*2-1.0
       ax(i) = z
    end do
    !!
    !! random angle in radian
    call random_number(z)
    z = z*2.0*3.141529
    !! generate rotation matrix
    call rmat_from_ax_angle(ax, z, rmat)
    !write(*,*) 'random rotation matrix:'
    !do i = 1, 3
     ! write(*,*) rmat(i,:)
    !end do
    do i = 1, nat
       coords(:,i) = matmul(rmat,coords(:,i))
    end do

    !!
    !! random translation
    !!
    !! random vector
    do i = 1, 3
       call random_number(z)
       z = z*2.0 - 1.0
       tr(i) = z
    end do
    !! random norm
    call random_number(z)
    z = z*10.0
    tr(:) = z*tr(:)/norm2(tr(:))
    ! write(*,*) 'random translate:'
    ! write(*,*) tr
    !! translate
    do i = 1, nat
       coords(:,i) = coords(:,i) + tr
    end do

    !!
    !! make reflection randomly
    !!
    call random_number(z)
    ! write(*,*) 'mirror:', nint(z)
    if( nint(z) .eq. 1 ) then
       !! make mirror
       coords(3,:) = -coords(3,:)
    endif

  end subroutine randomize

  subroutine rmat_from_ax_angle(ax_in, angle, rmat)
    !! wikipedia/rotation_matrix#quaternion
    implicit none
    real, dimension(3), intent(in) :: ax_in
    real, intent(in) :: angle
    real, dimension(3,3), intent(out) :: rmat

    real, dimension(3) :: ax
    real :: c, s, c_one, dr

    ax = ax_in/norm2(ax_in)

    c = cos(angle)
    s = sin(angle)
    c_one = 1.0 - c

    rmat(1,1) = ax(1)**2*c_one + c
    rmat(1,2) = ax(1)*ax(2)*c_one - ax(3)*s
    rmat(1,3) = ax(1)*ax(3)*c_one + ax(2)*s

    rmat(2,1) = ax(2)*ax(1)*c_one + ax(3)*s
    rmat(2,2) = ax(2)**2*c_one + c
    rmat(2,3) = ax(2)*ax(3)*c_one - ax(1)*s

    rmat(3,1) = ax(3)*ax(1)*c_one - ax(2)*s
    rmat(3,2) = ax(3)*ax(2)*c_one + ax(1)*s
    rmat(3,3) = ax(3)**2*c_one + c

    !call determinant(rmat,dr)
    !write(*,*) 'det r:',dr
  end subroutine rmat_from_ax_angle

  subroutine set_random_seed()
    ! ----- setting a random seed, based on current time -----
    integer :: i_seed
    integer, dimension(:), allocatable :: a_seed
    integer, dimension(1:8) :: dt_seed

    ! set random seed
    call random_seed(size=i_seed)
    allocate(a_seed(1:i_seed))
    call random_seed(get=a_seed)
    call date_and_time(values=dt_seed)
!    a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    a_seed(:)=dt_seed(8); a_seed(1)=dt_seed(8)
    call random_seed(put=a_seed)
    deallocate(a_seed)
  end subroutine set_random_seed

  subroutine determinant( a, d )
    !! determinant of a 3x3 matrix
    implicit none
    real, dimension(3,3), intent(in) :: a
    real, intent(out) :: d

    d = a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) + &
        a(1,2)*a(2,3)*a(3,1) - a(2,2)*a(3,1)*a(1,3) + &
        a(2,1)*a(3,2)*a(1,3) - a(3,3)*a(1,2)*a(2,1)

  end subroutine determinant


  subroutine random_permutation( n, list )
    !! generate a list of random indices of size n, such that
    !! no index repeats
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(out) :: list

    integer :: i, idx
    real :: z
    logical :: old

    !! initial values
    list(:) = 0

    do i = 1, n
       old = .true.
       do while( old )
          call random_number(z)
          !! generate index randomly in the range [1:n]
          idx = int( z*n ) + 1
          !! if this index already in list, skip and geenrate new
          if( any(list .eq. idx ) ) cycle
          !! if not, add it to list and stop loop for current index
          list(i) = idx
          old = .false.
       end do
    end do
  end subroutine random_permutation

