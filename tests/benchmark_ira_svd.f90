program benchmark_ira_svd
  use ira_mod, only: ira_rp, ira_svd
  implicit none

  character(len=512) :: path, arg, label
  integer :: unit, nat, niter, i, iter, ierr
  integer, allocatable :: typ1(:), typ2(:), perm(:), source_perm(:)
  real(ira_rp), allocatable :: coords1(:,:), coords2(:,:)
  real(ira_rp) :: angle, c, s, hd, rmsd, checksum
  real(ira_rp) :: rotation(3,3), generated_rotation(3,3), translation(3)
  integer :: clock_start, clock_stop, clock_rate

  call get_command_argument(1, path)
  if (len_trim(path) == 0) then
    error stop "usage: benchmark_ira_svd /path/to/100.xyz [iterations]"
  end if
  niter = 10000
  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) read(arg, *) niter

  open(newunit=unit, file=trim(path), status="old", action="read")
  read(unit, *) nat
  read(unit, *)
  allocate(typ1(nat), typ2(nat), perm(nat), source_perm(nat))
  allocate(coords1(3,nat), coords2(3,nat))
  do i = 1, nat
    read(unit, *) label, coords1(:,i)
    typ1(i) = 1
  end do
  close(unit)

  angle = 0.731_ira_rp
  c = cos(angle)
  s = sin(angle)
  generated_rotation = reshape([c, s, 0.0_ira_rp, &
                               -s, c, 0.0_ira_rp, &
                                0.0_ira_rp, 0.0_ira_rp, -1.0_ira_rp], [3,3])
  translation = [0.37_ira_rp, -0.91_ira_rp, 1.23_ira_rp]
  do i = 1, nat
    source_perm(i) = nat + 1 - i
    coords2(:,i) = matmul(generated_rotation, coords1(:,source_perm(i))) &
                   + translation
  end do
  typ2 = typ1(source_perm)

  checksum = 0.0_ira_rp
  call system_clock(clock_start, clock_rate)
  do iter = 1, niter
    call ira_svd(nat, typ1, coords1, nat, typ2, coords2, 1.2_ira_rp, &
                 rotation, translation, perm, hd, rmsd, ierr)
    if (ierr /= 0) error stop "ira_svd failed"
    checksum = checksum + hd + rmsd + rotation(1,1)
  end do
  call system_clock(clock_stop)

  write(*,'(a,i0)') "iterations=", niter
  write(*,'(a,f0.9)') "seconds=", &
    real(clock_stop-clock_start, ira_rp)/real(clock_rate, ira_rp)
  write(*,'(a,es24.16)') "checksum=", checksum
end program benchmark_ira_svd
