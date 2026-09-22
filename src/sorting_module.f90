module sorting_module

  !!
  !! Merge sort algorithm, obtained from Rosetta Code wiki on 29.11.2019:
  !!
  !!     https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort
  !!
  !! Slightly modified to sort 2D arrays in given input axis.
  !!
  use ira_precision
  implicit none

  private
  public :: sort
contains


  subroutine merge_a(A, B, C, ax)
    implicit none
    ! The targe attribute is necessary, because A .or. B might overlap with C.
    real(rp), target, intent(in) :: A(:,:), B(:,:)
    real(rp), target, intent(inout) :: C(:,:)
    integer(ip),      intent(in) :: ax

    integer(ip) :: i, j, k, sizea, sizeb, sizec

    sizea = size(A,2)
    sizeb = size(B,2)
    sizec = size(C,2)

    if (sizeA + sizeB > sizeC) stop "(1)"

    i = 1; j = 1
    do k = 1, sizeC
      if (i <= sizeA .and. j <= sizeB) then
        if (A(ax,i) <= B(ax,j)) then
          C(:,k) = A(:,i)
          i = i + 1
        else
          C(:,k) = B(:,j)
          j = j + 1
        end if
      else if (i <= sizeA) then
        C(:,k) = A(:,i)
        i = i + 1
      else if (j <= sizeB) then
        C(:,k) = B(:,j)
        j = j + 1
      end if
    end do
  end subroutine merge_a


  subroutine swap(ndim, x, y)
    implicit none
    integer(ip),               intent(in) :: ndim
    real(rp), dimension(ndim), intent(inout) :: x, y

    real(rp), dimension(ndim) :: tmp

    tmp = x
    x = y
    y = tmp
  end subroutine swap


  recursive subroutine mergesort(A, work, ax)
    implicit none
    real(rp),    intent(inout) :: A(:,:)
    real(rp),    intent(inout) :: work(:,:)
    integer(ip), intent(in) :: ax

    integer(ip) :: half, sizeA, ndim

    ndim = size(A,1)
    sizeA =size(A,2)
    half = (sizeA + 1) / 2
    if (sizeA < 2) then
      continue
    else if (sizeA == 2) then
      if (A(ax,1) > A(ax,2)) then
        call swap(ndim,A(:,1), A(:,2))
      end if
    else
      call mergesort(A(:, : half), work,ax)
      call mergesort(A(:, half + 1 :), work,ax)
      if (A(ax,half) > A(ax,half + 1)) then
        work(:,1 : half) = A(:,1 : half)
        call merge_a(work(:,1 : half), A(:,half + 1:), A, ax)
     end if
    end if
  end subroutine mergesort

  !> @brief Use mergesort to sort the input array by the chosen axis.
  subroutine sort(n, ndim, array, axis)
    implicit none
    integer(ip), intent(in) :: n
    integer(ip), intent(in) :: ndim
    real(rp), dimension(ndim, n), intent(inout) :: array
    integer(ip), intent(in) :: axis

    real(rp), allocatable :: work(:, :)

    !! allocate working space for the sorting routine
    allocate (work(1:ndim, 1:(n + 1)/2), source=0.0_rp)

    call mergesort(array, work, axis)

    deallocate (work)

  end subroutine sort

end module sorting_module
