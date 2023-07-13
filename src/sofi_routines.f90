!!
!! Copyright (C) 2023, MAMMASMIAS Consortium
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

  subroutine sofi_compute_all( nat, typ, coords, sym_thr, &
                               nmat, mat_list, perm_list, &
                               op_list, n_list, p_list, &
                               ax_list, angle_list, dmax_list, pg )
    !!
    !! compute all data from SOFI in single routine.
    !! coords are assumed to be already shifted to desired origin on input
    !!
    !! NOTE: all arrays have size nmax, but the actual output values go only up to index nmat,
    !! values beyond that index can be random.
    !!
    !!===== input variables
    !! nat                        :: number of atoms
    !! typ, dimension ( nat )     :: integer atomic types
    !! coords, dimension( 3, nat) :: positions of atoms
    !! sym_thr                    :: threshold for finding symmetries
    !!                               (not taken into account when making combinations)
    !!===== output variables
    !! nmat       :: number of symmetries found
    !! mat_list   :: symmetry matrices
    !! perm_list  :: permutation of atoms after applying each symmetry
    !! op_list    :: Character "Op" from the Schoenflies notation: Op n^p
    !!                (Id = identity, I = inversion, C = rotation, S = (roto-)reflection )
    !! n_list     :: Schoenflies n value
    !! p_list     :: Schoenflies p value
    !! ax_list    :: axis of operation of each symmetry
    !! angle_list :: angle of each symmetry, in units of 1/2pi,
    !!                i.e. angle=0.333 is 1/3 of full circle, or 120 degrees
    !! dmax_list  :: max difference of atomic positions of before/after symm transformation
    !! pg         :: name of Point group, e.g. D6h
    !!
    use sofi_tools, only: nmax
    implicit none
    !! ===== input
    integer,                 intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ
    real, dimension(3,nat),  intent(in) :: coords
    real,                    intent(in) :: sym_thr
    !! ===== output
    integer,                       intent(out) :: nmat
    real, dimension(3,3,nmax),     intent(out) :: mat_list
    integer, dimension(nat, nmax), intent(out) :: perm_list
    character(len=2), dimension(nmax), intent(out) :: op_list
    integer, dimension(nmax),      intent(out) :: n_list
    integer, dimension(nmax),      intent(out) :: p_list
    real, dimension(3, nmax),      intent(out) :: ax_list
    real, dimension(nmax),         intent(out) :: angle_list
    real, dimension(nmax),         intent(out) :: dmax_list
    character(len=10),             intent(out) :: pg

    real :: dum
    integer :: i
    real, dimension(3) :: rdum
    real, dimension(3,3) :: rmat
    logical :: verb

    !! get symmetries with sym_thr
    call sofi_get_symmops( nat, typ, coords, sym_thr, nmat, mat_list )

    !! get combos, can produce symm above sym_thr
    call sofi_get_combos( nat, typ, coords, nmat, mat_list )

    !! get ops, and unique angles and axes
    call sofi_unique_ax_angle( nmat, mat_list, op_list, ax_list, angle_list )

    !! get permuations and dmax
    call sofi_get_perm( nat, typ, coords, nmat, mat_list, perm_list, dmax_list )

    !! get name of pg
    verb = .false.
    call sofi_get_pg( nmat, mat_list, pg, verb )

    !! get the p, n values
    do i = 1, nmat
       rmat = mat_list(:,:,i)
       call sofi_analmat( rmat, op_list(i), n_list(i), p_list(i), rdum, dum )
    end do

  end subroutine sofi_compute_all


  subroutine sofi_struc_pg( nat, typ_in, coords_in, sym_thr, pg, verb )
    !! wrapper routine for directly getting only the PG name from a structure
    !! The structure is shifted to geometrical center!
    use sofi_tools, only: nmax
    implicit none
    integer,                 intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ_in
    real, dimension(3,nat),  intent(in) :: coords_in
    real,                    intent(in) :: sym_thr
    character(len=10),       intent(out) :: pg
    logical,                 intent(in) :: verb
    ! local
    integer :: n_op
    integer, allocatable :: perm_list(:,:)
    real, allocatable :: op_list(:,:,:)
    real, dimension(3) :: gc
    integer :: i
    integer, dimension(nat) :: typ
    real, dimension(3,nat) :: coords

    allocate( op_list(1:3,1:3,1:nmax) )
    allocate( perm_list(1:nat, 1:nmax) )

    !! working copies
    typ = typ_in
    coords = coords_in

    !! recenter to gc
    gc = sum( coords(:,:),2)/nat
    do i = 1, nat
       coords(:,i) = coords(:,i) - gc
    end do

    !! get list of symm operations
    call sofi_get_symmops( nat, typ, coords, sym_thr, n_op, op_list )
    ! call sofi_get_symmops( nat, typ, coords, sym_thr, n_op, op_list, perm_list )

    !! get PG name
    ! call sofi_get_pg( n_op, op_list, pg, verb )

    !! combos with higher sym_thr: find the missing Ops in case of '+' or '-' PGs
    call sofi_get_combos( nat, typ, coords, n_op, op_list )
    ! call sofi_get_combos( nat, typ, coords, n_op, op_list, perm_list )

    !! get new PG name
    call sofi_get_pg( n_op, op_list, pg, verb )

    deallocate( op_list, perm_list )

  end subroutine sofi_struc_pg


  subroutine sofi_get_symmops( nat, typ_in, coords_in, sym_thr, n_so, op_list )
  ! subroutine sofi_get_symmops( nat, typ_in, coords_in, sym_thr, n_so, op_list, perm_list )
  !! main SOFI routine
    !!====================================
    !!
    !! - operations from op_list can be applied for example the N-th entry by matmul():
    !!
    !!    theta(3,3) = op_list(:, :, N)
    !!    do i = 1, natoms
    !!       coords_rot(:,i) = matmul( theta, coords_orig(:,i) )
    !!    end do
    !!
    !! - and adding permutations from perm_list:
    !!
    !!    coords_rot(:, perm_list(:, N) ) = coords_rot(:,:)
    !!
    !! - after these 2 operations, coords_rot and coords_orig should be equal
    !!   atom-by-atom ( or within sym_thr )
    !!
    use sofi_tools, only: nmax, m_thr
    implicit none
    integer,                 intent(in) :: nat        ! number of atoms
    integer, dimension(nat), intent(in) :: typ_in     ! atomic types
    real, dimension(3,nat),  intent(in) :: coords_in  ! atomic xyz
    real,                    intent(in) :: sym_thr    ! threshold for symmetry (in terms of hausdorff dist.)
    !!
    integer,                           intent(out) :: n_so       ! number of found symm operations
    real, dimension(3,3,nmax),         intent(out) :: op_list   ! the list of operations
    ! integer, dimension(1:nat, 1:nmax), intent(out) :: perm_list  ! list of associated permutations

    !! local
    integer, dimension(nat) :: typ
    real, dimension(3,nat) :: coords
    real, allocatable :: d_o(:,:)
    integer :: i, j, k, l, m, mm
    real, dimension(3,3) :: theta, beta, gamma
    ! real, dimension(3) :: ax, ax1
    logical :: fail1
    real :: dd, d_i, d_j, dh
    integer :: nbas, nn
    integer :: ti, tj
    real :: small_norm
    ! real, dimension(3) :: v1, v2

    ! if( size(op_list,3) .ne. nmax .or. size(perm_list,2) .ne. nmax ) then
    if( size(op_list,3) .ne. nmax ) then
    write(*,*) 'the op_list on input should have maxsize:',nmax
       stop
    end if

    !! zero the output data
    nbas = 0
    op_list(:,:,:) = 0.0
    ! perm_list(:,:) = 0

    !! set local copy
    typ(:) = typ_in(:)
    coords(:,:) = coords_in(:,:)

    !! sort by d
    allocate( d_o(1:2,1:nat))
    do i = 1, nat
      d_o(1,i) = norm2(coords(:,i))
      d_o(2,i) = real(i)
    end do

    !! sorting routines from IRA lib
    call sort( nat, 2, d_o, 1)
    coords = coords(:,nint(d_o(2,:)))
    typ = typ( nint(d_o(2,:)))

    !! find smallest atom-atom dist, for first cshda thr
    !!=============
    !! The idea is to be a bit "flexible" in first cshda, to
    !! allow some margin of error which hopefully disappears in refine.
    !! This could be played a bit and optimized, i think max dd could
    !! be proportional to half the smallest atom-atom distance.
    dd = 99.0
    do i = 1, nat
      do j = i+1, nat
          dd = min( dd, norm2(coords(:,i) - coords(:,j)))
      end do
    end do


    !! add id
    theta(:,:) = 0.0
    do i = 1, 3
      theta(i,i) = 1.0
    end do

    call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, dh, 0.5 )
    ! call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, perm_list, dh, 0.5 )

    !! try inversion
    theta(:,:) = 0.0
    do i = 1, 3
      theta(i,i) = -1.0
    end do

    call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, dh, 0.5 )
    ! call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, perm_list, dh, 0.5 )

    ! do i = 1, nat
    !    write(*,*) i, norm2( coords(:,i))
    ! end do


    !! set beta
    small_norm = 1e-1
    small_norm = max(sym_thr, 1e-1)
    ! small_norm = sym_thr
    ! small_norm = 5e-1
    do i = 1, nat
      if( norm2( coords(:,i)) .lt. small_norm) cycle
      do j = 1, nat
          if( norm2( coords(:,j)) .lt. small_norm) cycle
          !!
          !! call to IRA lib
          call set_orthonorm_bas( coords(:,i), coords(:,j), beta, fail1 )
          if( fail1 ) cycle
          !!
          exit
      end do
      exit
    end do

    ! write(*,*) i,j
    if( i .gt. nat .or. j .gt. nat ) then
      write(*,*) repeat('%',40)
      write(*,*) "ERROR: cannot set beta. Linear structure?"
      write(*,*) repeat('%',40)
      stop
    endif

    !! distances of atoms i,j
    d_i = norm2( coords(:,i))
    d_j = norm2( coords(:,j))
    ti = typ(i)
    tj = typ(j)
    ! write(*,*) 'dd',dd
    ! write(*,*) d_i, d_j, i, j

    ! write(*,*) 'beta is'
    ! write(*,'(3f9.4)') beta(1,:)
    ! write(*,'(3f9.4)') beta(2,:)
    ! write(*,'(3f9.4)') beta(3,:)

    ! write(*,*) nat
    ! write(*,*) "properties=species:S:1:pos:R:3:id:I:1"
    ! do i = 1, nat
    !    write(*,*) typ(i), matmul((beta), coords(:,i)), i
    !    ! write(*,*) typ(i), coords(:,i), i
    ! end do


    !!========
    !! main loop
    !!========
    nn = 0
    do k = 1, nat
       if( abs( d_i - norm2(coords(:,k))) .gt. 1.1*sym_thr ) cycle
       if( typ(k) .ne. ti ) cycle
       do l = 1, nat
          if( k.eq.l) cycle
          if( abs( d_j - norm2(coords(:,l))) .gt. 1.1*sym_thr ) cycle
          if( typ(l) .ne. tj ) cycle
          !!
          !! call to IRA lib
          !! (re)-set gamma
          call set_orthonorm_bas( coords(:,k), coords(:,l), gamma, fail1)
          if( fail1 ) cycle
          !!
          do mm = 1, 5
            m = mm - 4
            if( m .eq. 0 ) cycle
            nn = nn + 1
            !!
            !! call to IRA lib
            !! set gamma
            call set_orthonorm_bas( coords(:,k), coords(:,l), gamma, fail1)
            !! flip ax of gamma
            gamma(abs(m),:) = sign(1,m)*gamma(abs(m),:)
            !! set trial theta
            theta = matmul( transpose(beta), gamma)
            !!
            ! write(*,*) repeat('-',20)
            ! write(*,*) "generated next mat"
            ! write(*,*) theta(1,:)
            ! write(*,*) theta(2,:)
            ! write(*,*) theta(3,:)
            ! call sofi_analmat( theta, ... )
            call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, dh, m_thr )
            ! call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, perm_list, dh, m_thr )
            ! write(*,*) 'op exiting try_sofi:'
            ! call sofi_analmat( theta, ... )
            ! write(*,'(2i4,3f9.5,x,f9.4,a5,f9.4)') k,l,ax1,d1,'dh:',dh
            ! theta = matmul(beta,transpose(gamma))
            ! call try_sofi( theta, nat, typ, coords, sym_thr, dd, nbas, op_list, perm_list, dh, m_thr )
          end do
       end do
    end do
    ! write(*,*)
    ! write(*,*) 'n tests, nn =',nn
    ! write(*,*) 'loop 0',nbas

    !! For the combinations of SymmOps we should
    !! lower the matrix distance threshold: value 0.35 captures 1/24*2pi rotations,
    !! this should be good-enough for regular PGs. The groups C, S, D, T, O, I have
    !! the smallest separation of elements at 1/12*2pi, which is at m_thr=~0.73
    !!==========
    !! NOTE: the above idea is nonsense, since combinations cannot produce smaller
    !! angle operation than original..... Therefore, the original m_thr needs to be low!
    !!
    ! m_thr = 0.35
    ! m_thr = 0.5
    !! combos
    ! do m = 1, 10
    !   ii = nbas
    !   do i = 1, nbas
    !       do j = 1, nbas
    !         ! write(*,*) 'combo',i,j
    !         theta = matmul( op_list(:,:,i), op_list(:,:,j))
    !         call try_sofi( theta, nat, typ, coords, sym_thr, dd, ii, op_list, perm_list, dh, m_thr )
    !         ! write(*,*) i,j,dh
    !         ! write(*,'(3f9.4)') theta(1,:)
    !         ! write(*,'(3f9.4)') theta(2,:)
    !         ! write(*,'(3f9.4)') theta(3,:)
    !       end do
    !   end do
    !   write(*,*) 'loop',m,ii, nbas
    !   if( nbas .eq. ii ) exit
    !   nbas = ii
    ! end do

    !! set output number of operations
    n_so = nbas

    ! write(*,*) 'exiting get_symmops'
    deallocate( d_o )

  end subroutine sofi_get_symmops


  subroutine sofi_get_perm( nat, typ, coords, nbas, bas_list, perm_list, dmax_list )
    implicit none
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ
    real, dimension(3,nat), intent(in) :: coords
    integer, intent(in) :: nbas
    real, dimension(3,3,nbas), intent(in) :: bas_list
    integer, dimension(nat, nbas), intent(out) :: perm_list
    real, dimension(nbas), intent(out) :: dmax_list

    integer :: i, j
    real, dimension(3,3) :: rmat
    integer, dimension(nat) :: t_local
    real, dimension(3,nat) :: c_local
    integer, dimension(nat) :: found
    real, dimension(nat) :: dists

    do i = 1, nbas
       !! get matrix from list
       rmat = bas_list(:,:,i)
       !! transform struc
       do j = 1, nat
          t_local(j) = typ(j)
          c_local(:,j) = matmul(rmat, coords(:,j) )
       end do
       !! compute cshda with high thr (accept anything)
       call cshda( nat, typ, coords, nat, t_local, c_local, 99.9, found, dists )
       ! if( any(found) .eq. 0 ) then
       !    write(*,*) "!! ERROR IN: sofi_get_perm, failed cshda?"
       !    stop
       ! endif
       perm_list(:,i) = found(1:nat)
       dmax_list(i) = maxval( dists(1:nat) )
       ! write(*,*) i, maxval(dists(1:nat))
       ! write(*,*) rmat(1,:)
       ! write(*,*) rmat(2,:)
       ! write(*,*) rmat(3,:)

    end do

  end subroutine sofi_get_perm


  subroutine sofi_get_combos( nat, typ, coords, nbas, bas_list )
  ! subroutine sofi_get_combos( nat, typ, coords, nbas, bas_list, perm_list )
  !! find all combos with high sym_thr: this will accept anything.
    !! However, some strange matrices cannot be constructed from combos of
    !! matrices which are already found as symmetries with sym_thr value.
    use sofi_tools, only: m_thr, nmax
    implicit none
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: typ
    real, dimension(3, nat), intent(in) :: coords
    integer, intent(inout) :: nbas
    real, dimension(3, 3, nmax), intent(inout) :: bas_list
    ! integer, dimension(nat, nmax), intent(inout) :: perm_list

    integer :: m, i, j, ii
    real, dimension(3,3) :: theta
    real :: dh, dd, sym_thr

    !! hard-coded thresholds (accept anything)
    dd = 10.0
    sym_thr = 5.0

    ! write(*,*) 'in sofi_get_combos'
    ! do i = 1, nbas
    !    write(*,*) i
    !    write(*,*) bas_list(1,:,i)
    !    write(*,*) bas_list(2,:,i)
    !    write(*,*) bas_list(3,:,i)
    ! end do

    do m = 1, 10
       ii = nbas
       do i = 1, nbas
          do j = 1, nbas
             ! write(*,*) 'combo',i,j
             theta = matmul( bas_list(:,:,i), bas_list(:,:,j))
             call try_sofi( theta, nat, typ, coords, sym_thr, dd, ii, bas_list, dh, m_thr )
             ! call try_sofi( theta, nat, typ, coords, sym_thr, dd, ii, bas_list, perm_list, dh, m_thr )
             ! if( dh .lt. sym_thr ) then
             !    write(*,*) i,j,dh
             !    write(*,'(3f9.4)') theta(1,:)
             !    write(*,'(3f9.4)') theta(2,:)
             !    write(*,'(3f9.4)') theta(3,:)
             ! endif

          end do
       end do
       if( nbas .eq. ii ) exit
       ! write(*,*) 'loop',m,ii, nbas
       nbas = ii
    end do

  end subroutine sofi_get_combos



  subroutine try_sofi( theta, nat, typ_in, coords_in, sym_thr, dd, nbas, op_list, dh, m_thr )
  ! subroutine try_sofi( theta, nat, typ_in, coords_in, sym_thr, dd, nbas, op_list, perm_list, dh, m_thr )
  !! try if the matrix theta gives dH within 5*sym_thr, then send it to refine.
    !! If matrix is valid and new, add it to op_list, and perm_list
    use sofi_tools, only: nmax
    implicit none
    real, dimension(3,3),          intent(inout) :: theta
    integer,                       intent(in) :: nat
    integer, dimension(nat),       intent(in) :: typ_in
    real, dimension(3,nat),        intent(in) :: coords_in
    real,                          intent(in) :: sym_thr
    real,                          intent(in) :: dd !! thr for first cshda, "smallest atom-atom distance"
    integer,                       intent(inout) :: nbas
    real, dimension(3,3,nmax),     intent(inout) :: op_list
    ! integer, dimension(nat, nmax), intent(inout) :: perm_list
    real,                          intent(out) :: dh
    real,                          intent(in) ::m_thr

    integer, dimension(nat) :: typ
    real, dimension(3,nat) :: coords
    integer :: i
    real, allocatable :: dists(:)
    integer, allocatable :: found(:), perm(:)
    logical :: not_crazy, is_new
    logical :: is_valid

    ! write(*,*) '>> entering try_sofi'
    !! local copy
    typ = typ_in
    coords = coords_in

    dh = 99.9
    is_valid = .false.

    !! apply theta
    do i = 1, nat
      coords(:,i) = matmul( theta, coords(:,i) )
    end do

    !! check if theta is already known
    call is_new_sofi( theta, nbas, op_list, m_thr, is_new )
    ! write(*,*) 'is new',is_new
    if( .not. is_new ) return

    !! call to IRA lib
    !! get first cshda
    allocate( found(1:nat))
    allocate( dists(1:nat))
    allocate( perm(1:nat))
    call cshda( nat, typ_in, coords_in, &
        nat, typ, coords, &
        1.1*dd, found, dists )
    dh = maxval(dists,1)
    ! write(*,'(x,a,x,f9.4)') 'initial dh',dh
    if( dh .gt. 5.0*sym_thr) return
    ! do i = 1, nat
    !    write(*,*) i, found(i), dists(i)
    ! end do
    ! write(*,*) '1st dh:',dh, sum(dists**2)



    not_crazy = .true.
    !! cshda returned no assignments
    if( any(found .eq. 0)) not_crazy = .false.

    ! write(*,*) 'not_crazy',not_crazy
    if( .not. not_crazy ) return
    perm = found

    !! permute
    typ = typ(found)
    coords(:,:) = coords(:,found)
    ! write(*,*) nat
    ! write(*,*) 'ref'
    ! do i = 1, nat
    !    write(*,*) typ_in(i), coords_in(:,i)
    ! end do
    ! write(*,*) nat
    ! write(*,*) 'me'
    ! do i = 1, nat
    !    write(*,*) typ(i), coords(:,i)
    ! end do


    if( not_crazy) then

      !! refine
      call refine_sofi( nat, typ_in, coords_in, nat, typ, coords, theta, dh, perm )

      !! dh after refinement is good
      if( dh .le. sym_thr ) is_valid = .true.

    endif

    ! write(*,*) 'is valid:',is_valid, dh, sym_thr

    if( is_valid ) then
      !! check again if new (could change due to refine)
      call is_new_sofi( theta, nbas, op_list, m_thr, is_new )
      if( is_new ) then
          ! write(*,*) ':: adding'
          ! perm = found
          ! call add_sofi( nat, theta, perm, nbas, op_list, perm_list )
          call add_sofi( nat, theta, nbas, op_list )
      end if
    end if


    deallocate( found, dists )
    deallocate( perm )
    ! write(*,*) '::: exit try_sofi'
  end subroutine try_sofi


  subroutine refine_sofi( nat_ref, typ_ref, coords_ref, nat, typ_in, coords_in, theta, dh, perm )
    !! few steps of ICP-like algo...
    !! NOTE: assume atoms are already permuted at input
    !! *_ref is static, original structure
    !!
    implicit none
    integer,                     intent(in) :: nat_ref
    integer, dimension(nat_ref), intent(in) :: typ_ref
    real, dimension(3,nat_ref),  intent(in) :: coords_ref
    integer,                     intent(in) :: nat
    integer, dimension(nat),     intent(in) :: typ_in
    real, dimension(3,nat),      intent(in) :: coords_in
    real, dimension(3,3),        intent(inout) :: theta
    real,                        intent(inout) :: dh
    integer, dimension(nat),     intent(out) :: perm

    integer :: i, n, k
    integer, dimension(nat) :: found, found_new
    real, dimension(nat) :: dists
    integer, dimension(nat) :: typ
    real, dimension(3,nat) :: coords
    real, dimension(3,3) :: svd_r
    real, dimension(3) :: svd_t  !, ax
    ! real :: dd, det1

    !! working copies
    typ = typ_in
    coords = coords_in

    !! max steps
    n = 3

    ! write(*,*) 'theta on refine input:'
    ! write(*,'(3f9.4)') theta

    !! found is 1,2,3, ..
    found(:) = (/(i,i=1,nat)/)

    do i = 1, n

      ! write(*,*) i
      ! write(*,'(3f9.4)') theta

      !! call to IRA lib
      !! find correction from svd
      call svd_forcerot( nat_ref, typ_ref, coords_ref, &
            nat, typ, coords, svd_r, svd_t)

      !! this shoud not happen, positive determinant is forced by svd_forcerot
      ! call determinant3x3(svd_r,det1)
      ! if( det1 .lt. -0.5 ) svd_r(3,:) = -svd_r(3,:)


      !! only for ptinting, check what is this matrix
      ! call sofi_analmat( theta, ... )
      ! write(*,'(a,x,f12.6)') 'svd angle:',dd*360
      ! write(*,*) 'svd_r on step',i
      ! write(*,'(3f9.4)') svd_r

      ! write(*,*) 'svd_t'
      ! write(*,'(3f9.4)') svd_t

      !! apply
      do k = 1, nat
          coords(:,k) = matmul(svd_r,coords(:,k)) + svd_t
      end do
      ! write(*,*) nat
      ! write(*,*) 'me in svd, rotated'
      ! do j = 1, nat
      !    write(*,*) typ(j), coords(:,j)
      ! end do

      !! add correction to theta
      theta = matmul(svd_r, theta)

      !! call to IRA lib
      !! get permutation, with high thr (accept anything at this point)
      call cshda( nat_ref, typ_ref, coords_ref, &
            nat, typ, coords, &
            99.0, found_new, dists )

      ! write(*,*) 'found, dists'
      ! do k = 1, nat
      !    write(*,*) k, found(k), dists(k)
      ! end do
      ! write(*,*) 'dh',maxval(dists,1), sum(dists**2)

      ! write(*,*) 'cor:',i,maxval(dists)

      !! permute
      typ(:) = typ( found_new )
      coords(:,:) = coords(:, found_new )
      perm = perm(found_new)
      dh = maxval(dists,1)
      ! write(*,*) i, dh

      !! if no permutation is found, the refinement is finished
      if( all(found_new .eq. found) ) exit

    end do

    ! write(*,*) 'theta on refine output:'
    ! write(*,'(3f9.4)') theta


  end subroutine refine_sofi


  subroutine is_new_sofi( rmat, nbas, op_list, m_thr, is_new )
  !! check if rmat is new in op_list
    use sofi_tools, only: nmax, matrix_distance
    implicit none
    real, dimension(3,3),      intent(in) :: rmat
    integer,                   intent(in) :: nbas
    real, dimension(3,3,nmax), intent(in) :: op_list
    real,                      intent(in) :: m_thr
    logical,                   intent(out) :: is_new

    integer :: i
    logical :: eq_perm, eq_rmat
    real :: dd
    real :: det1, det2, ddet
    real :: matrix_thr

    !! Currently it is decided to compare rmat matrix to all the previously found
    !! matrices to decide whether rmat is new or not.
    !! Alternatively, this decision could be to compare permutation 'perm' to all the
    !! permutations in 'perm_list', since each symmOp should have a unique permutation.
    !! Or not? Mirror triangle over the plane gives same permutation as orig...


    !! threshold on equivalence of matrices, higher value can speed up overall calc,
    !! but too high can think different matrices are equal.
    !! Value around 0.5 corresponds to a
    !! rotation of matrix by cca 1/18*2pi, which seems generally ok.
    !! Should this be let as input?! Done.
    ! matrix_thr = 0.5
    matrix_thr = m_thr

    !! initial values
    is_new = .true.
    eq_perm = .false.
    eq_rmat = .false.

    do i = 1, nbas
      !! if rmat equal
      call matrix_distance( rmat, op_list(:,:,i), dd)
      if( dd .lt. matrix_thr ) then
          !! call to IRA
          !! check if determinants equal (should be, but maybe not when high matrix_thr)
          call determinant3x3( rmat, det1 )
          call determinant3x3( op_list(:,:,i), det2)
          ddet = abs( det1 - det2 )
          ! write(*,*) 'dmat:', dd
          ! write(*,*) 'det1:', det1, 'det2:',det2

          !! value of det can be +1 or -1, so no need for flexible thr here
          if( ddet .lt. 1e-3 ) then
            eq_rmat = .true.
          endif

          ! write(*,*) 'rmat'
          ! write(*,'(3f9.4)') rmat
          ! write(*,*) 'op_list(:,:,i)'
          ! write(*,'(3f9.4)') op_list(:,:,i)
      endif

      ! is_new = .not.(eq_perm .and. eq_rmat)
      is_new = .not.eq_rmat

      if( .not. is_new ) return
    end do

    ! write(*,*) "is new:", is_new
  end subroutine is_new_sofi


  subroutine add_sofi( nat, rmat, nbas, mat_list )
  ! subroutine add_sofi( nat, rmat, perm, nbas, op_list, perm_list )
  !! add rmat into mat_list, add perm in perm_list
    use sofi_tools, only: nmax
    implicit none
    integer,                      intent(in) :: nat
    real, dimension(3,3),         intent(in) :: rmat
    ! integer, dimension(nat),      intent(in) :: perm
    integer,                      intent(inout) :: nbas
    real, dimension(3,3,nmax),    intent(inout) :: mat_list
    ! integer, dimension(nat,nmax), intent(inout) :: perm_list


    !! increment nbas
    nbas = nbas + 1
    if( nbas .gt. nmax ) then
      write(*,*) "subroutine add_sofi::: ERROR ADDING RMAT, LIST TOO SMALL"
      return
    end if

    mat_list(:,:,nbas) = rmat
    ! perm_list(:,nbas) = perm

  end subroutine add_sofi

  subroutine sofi_get_pg( nbas, op_list, pg, verb )
    !! flowchart to determine PG notation from op_list
    !! online: https://symotter.org/assets/flowchart.pdf
    !!
    use sofi_tools
    implicit none
    integer,                   intent(in) :: nbas
    real, dimension(3,3,nbas), intent(in) :: op_list
    character(len=10),         intent(out) :: pg
    logical,                   intent(in) :: verb

    real, dimension(3) :: ax, cn_ax
    real, dimension(3,3) :: rmat
    real :: angle, dd, dot, cross
    integer :: i, j
    real :: pi
    character(len=2), dimension(nbas) :: op
    integer, dimension(nbas) :: n_int, power
    real, dimension(3,nbas) :: ax_list
    integer, dimension(nbas) :: skip_ax
    integer, dimension(nbas) :: multip_ax !! total number of operations on the axis of this op
    integer :: max_n_val, max_n_loc
    integer :: nax, nr_n
    real, dimension(5,nbas) :: uniq_ax  !! 4th dim is largest Cn on this ax, 5th is counter on op
    logical :: has_inversion, has_sigma, has_cn, has_sigma_h, has_sigma_v, has_s2n, c2_perp, has_sigma_d
    integer :: n_c2, n_in_plane, ndum
    real, dimension(3) :: cc
    integer :: k, pg_err
    real :: dotj, dotk

    if( verb ) then
       write(*,*) repeat('=',20)
       write(*,*) "number of SymmOps entering get_pg:",nbas
    endif

    pi=4.0*atan(1.0)

    !! start from below
    pg = 'C1'
    if( nbas .eq. 1 ) return

    has_inversion = .false.
    has_sigma = .false.
    has_cn = .false.
    has_sigma_h = .false.
    has_sigma_v = .false.
    has_sigma_d = .false.
    has_s2n = .false.
    c2_perp = .false.

    !! classify all op_list
    do i = 1, nbas
       rmat = op_list(:,:,i)
       call sofi_analmat( rmat, op(i), n_int(i), power(i), ax_list(:,i), angle )
       if(verb) write(*,'(i2,a3,i3,2x,3f8.3,x,f9.4)') i, op(i), n_int(i), ax_list(:,i), angle
    end do

    uniq_ax(:,:) = 0.0
    nax = 0
    multip_ax(:) = 0
    skip_ax(:) = 0


    !! some basics
    has_inversion = find_inversion( nbas, op )
    has_sigma = find_sigma( nbas, op, n_int )
    has_cn = find_cn( nbas, op, n_int )


    !! group the Ops of the list by axis
    do i = 1, nbas
       if( op(i) .eq. 'I' ) cycle
       if( op(i) .eq. 'Id' ) cycle
       !!
       if( skip_ax(i) .eq. 1 ) cycle
       ax = ax_list(:,i)
       if(verb) write(*,'(a4,3f9.4)') 'ax:',ax
       !! add to uniq_ax
       ! multip_ax(i) = multip_ax(i) + 1
       nax = nax + 1
       uniq_ax(1:3,nax) = ax

       max_n_val = 0
       do j = 1, nbas
          !! skip Id and I
          if( op(j) == 'Id' )cycle
          if( op(j) == 'I' )cycle
          !!
          if( abs(dot_product(ax, ax_list(:,j))) .gt. 0.999 ) then
             !! is same ax
             if(verb) write(*,'(a3,g0,a1,g0)') op(j), n_int(j),'^', power(j)
             skip_ax(j) = 1
             !! keep maximal n value of C operations
             if( op(j) == 'C') max_n_val = max( max_n_val, n_int(j))
             !! count how many ops have this ax
             uniq_ax(5,nax) = uniq_ax(5,nax) + 1.0
          end if
       end do
       uniq_ax(4,nax) = real(max_n_val)
    end do


    !! put multiplicity of ax to each op
    do i = 1, nbas
       !! find which ax from uniq list
       do j = 1, nbas
          if( abs(dot_product(ax_list(:,i), uniq_ax(1:3,j))) .gt. 0.99 ) then
             multip_ax(i) = nint(uniq_ax(5,j))
             exit
          endif
       end do
    end do

    ! write(*,'(24i3)') nint( uniq_ax(5,:) )
    ! write(*,'(24i3)') multip_ax(:)
    ! do i = 1, nax
    !    write(*,'(3f9.4,x,i0)') uniq_ax(1:3,i), nint(uniq_ax(4,i))
    ! end do

    !! find largest n of C ops, there are at least 2 ops on list
    max_n_val = 0
    max_n_loc = 2
    if( count(op .eq. 'C') .gt. 0 ) then
       max_n_val = maxval( n_int(:), mask=(op .eq. 'C'))
       max_n_loc = maxloc( n_int(:), dim=1, mask=(op .eq. 'C'))
    endif
    !!
    !! special case for D2: n of principal ax is equal to n of other C ax
    !!
    !! D2 has all axes multip 1, any of the C2 axes can be principal,
    !! but D2d has 2 C2 axes which are not prncipal, with multip=1, and principal C2 with multip=3
    if( maxval(multip_ax) .gt. 1 ) then
       !! there are axes with multip > 1
       if( multip_ax( max_n_loc) == 1 ) then
          !! we have chosen the wrong axis, choose again among axes with multip > 1
          max_n_val = maxval( n_int(:), mask=(op .eq. 'C' .and. multip_ax .gt. 1))
          max_n_loc = maxloc( n_int(:), dim=1, mask=(op .eq. 'C' .and. multip_ax .gt. 1))
       endif
    end if

    !! can happen if there are no C operations in PG, for example Cs. Choose the ax of second op
    if( max_n_loc .eq. 0 ) max_n_loc = 2


    if(verb) then
      write(*,*) 'largest n:', max_n_val, max_n_loc
      write(*,*) 'ax:',ax_list(:,max_n_loc)
    endif

    !! get how many C ax have this n
    nr_n = count( nint(uniq_ax(4,:)) .eq. max_n_val )
    ! write(*,*) 'nr_n', nr_n

    ! write(*,*) 'has inversion:',has_inversion
    ! write(*,*) 'has sigma:', has_sigma
    ! write(*,*) 'has cn:', has_cn
    ! write(*,*) 'max_n_val:',max_n_val
    ! write(*,'(a14,3f9.4)') 'principal ax:',ax_list(:,max_n_loc)

    !! select ax with largest n
    cn_ax = ax_list(:,max_n_loc)


    !!
    !! flowchart, top part
    !!
    if( max_n_val .ge. 3 .and. nr_n .ge. 2 ) then
       !!
       !! is some T# group (n < 4)
       if( max_n_val .ge. 4) then
          if( max_n_val .ge. 5 ) then
             !! is I#
             pg = 'I'
             if( has_sigma .or. has_inversion) pg = 'Ih'
          else
             !! is O#
             pg = 'O'
             if( has_sigma .or. has_inversion) pg = 'Oh'
          endif
       else
          if( has_sigma ) then
             pg = 'Td'
             if( has_inversion ) pg = 'Th'
          else
             pg = 'T'
          endif
          !!
       endif
       !!
       ! return
       goto 100
    endif

    !!
    !! bottom part of flowchart
    !!
    !! see if n C2 ax are perp to cn_ax
    dd = 0.0
    n_c2 = 0
    ! write(*,*) '>> checking c2 perp'
    do i = 1, nbas
       if( op(i) .ne. 'C') cycle
       if( n_int(i) .ne. 2) cycle
       if( i .eq. max_n_loc ) cycle
       ax = ax_list(:,i)
       !! keep dot prod of all
       dot = abs( dot_product( cn_ax, ax))
       !! can be same ax as cn_ax
       if( dot .gt. 0.99) cycle
       ! write(*,*) op(i), n_int(i)
       ! write(*,'(5x,i2,f5.2,3f9.4)') i, dot, ax
       dd = dd + dot
       !! count c2 axes
       n_c2 = n_c2 + 1
       ! write(*,*) i, dd
    end do
    if( n_c2 .gt. 0) dd = dd/n_c2
    ! write(*,*) 'n_c2',n_c2
    ! write(*,*) 'dd is:',dd
    if( dd .lt. 0.01 .and. n_c2 .gt. 0 ) c2_perp = .true.
    ! write(*,*) '>> c2 are perp',c2_perp

    !! check if has sigma_h
    !! horizontal refl plane: axis of plane is parallel to principal, dot( ax_s0, principal ) = 1.0
    ! write(*,*) '>> checking sigma_h'
    do i = 1, nbas
       if( op(i) .ne. 'S') cycle
       if( n_int(i) .ne. 0) cycle
       dot = abs( dot_product(cn_ax, ax_list(:,i)))
       if( dot .gt. 0.99 ) has_sigma_h = .true.
       ! write(*,*) i, dot
    end do
    ! write(*,*) '>> sigma_h found:',has_sigma_h



    !! check if has n sigma_v:
    !! vertical plane, contains the principal ax: the norm of cross( ax_s0, principal ) is 1.0
    ! write(*,*) '>> checking sigma_v'
    n_in_plane = 0
    do i = 1, nbas
       if( op(i) .ne. 'S') cycle
       if( n_int(i) .ne. 0) cycle
       !! cross prod( ax, cn_ax ) = 1 for ax in plane
       call cross_prod( ax_list(:,i), cn_ax, cc)
       cross=norm2(cc)
       ! write(*,'(5x,i2,f5.2,3f9.4)') i, cross, ax_list(:,i)
       if( abs( cross ) .gt. 0.99 ) n_in_plane = n_in_plane + 1
    end do
    ! write(*,*) 'n_in_plane',n_in_plane, n_c2
    ! if( n_in_plane .eq. n_c2) has_sigma_v = .true.
    if( n_in_plane .gt. 0) has_sigma_v = .true.
    ! write(*,*) '>> sigma_v found:',has_sigma_v


    !! check if has sigma_d
    !! diagonal plane: contains principal axis, and bisects
    !! the angle between 2 C2 axes that are perp to principal.

    !! Condition checked here: find s0 where dot( ax_s0, sum(ax_C2i, ax_C2j) ) = 0
    !! This condition is true for the "absolute direction" of sum( ax_C2i, ax_C2j),
    !! in the sense that both +/- directions need to be checked for both i,j

    ndum = 0
    ! write(*,*) '>> checking sigma_d'
    !! "contains principal axis" is equal to has_sigma_v:
    if( has_sigma_v ) then

       !! first, find two C2 axes perpendicular to principal ax: i,j
       lp:do i = 1, nbas
          if( op(i) .ne. 'C' ) cycle
          if( n_int(i) .ne. 2) cycle

          dot = abs( dot_product( ax_list(:,i), cn_ax))
          !! i is perp to cn
          if( dot .lt. 0.01) then
             !! select another
             do j = 1, nbas
                if( j .eq. i ) cycle
                if( op(j) .ne. 'C' ) cycle
                if( n_int(j) .ne. 2) cycle

                dotj = abs( dot_product( ax_list(:,j), cn_ax))
                if( dotj .lt. 0.01) then
                   ! write(*,*) i, j, dot,dotj
                   ! write(*,*) ax_list(:,i)
                   ! write(*,*) ax_list(:,j)
                   !! have 2 ax perp to cn, make sum, and dot with s0.
                   !! need to do this for all +/- combinations of i,j,
                   !! because their ax can have any orientation +/-
                   !! + +
                   ax = ax_list(:,i) + ax_list(:,j)
                   ax = ax/norm2(ax)
                   ! write(*,*) ax
                   !! dot with s0 should be zero -> s0 bisects axes i,j
                   do k = 1, nbas
                      !! loop through s0 ops
                      if( op(k) .ne. 'S') cycle
                      if( n_int(k) .ne. 0) cycle
                      !! dot with ax_s0
                      dotk = abs( dot_product(ax, ax_list(:,k)) )
                      ! write(*,*) i,j,k, dotk
                      !! check if this s0 plane contains principal
                      call cross_prod( ax_list(:,k), cn_ax, cc)
                      cross=norm2(cc)
                      if( dotk .lt. 0.01 .and. abs(cross) .gt. 0.99 ) then
                         has_sigma_d = .true.
                         exit lp
                      endif
                   end do
                   !! + -
                   ax = ax_list(:,i) - ax_list(:,j)
                   ax = ax/norm2(ax)
                   ! write(*,*) ax
                   !! dot with s0 should be zero -> s0 bisects axes i,j
                   do k = 1, nbas
                      if( op(k) .ne. 'S') cycle
                      if( n_int(k) .ne. 0) cycle
                      dotk = abs( dot_product(ax, ax_list(:,k)) )
                      call cross_prod( ax_list(:,k), cn_ax, cc)
                      cross=norm2(cc)
                     ! write(*,*) i,j,k, dotk
                      if( dotk .lt. 0.01 .and. cross.gt.0.99 ) then
                         has_sigma_d = .true.
                         exit lp
                      endif
                   end do
                   !! - +
                   ax = -ax_list(:,i) + ax_list(:,j)
                   ax = ax/norm2(ax)
                   ! write(*,*) ax
                   !! dot with s0 should be zero -> s0 bisects axes i,j
                   do k = 1, nbas
                      if( op(k) .ne. 'S') cycle
                      if( n_int(k) .ne. 0) cycle
                      dotk = abs( dot_product(ax, ax_list(:,k)) )
                      call cross_prod( ax_list(:,k), cn_ax, cc)
                      cross=norm2(cc)
                      ! write(*,*) i,j,k, dotk
                      if( dotk .lt. 0.01 .and. cross.gt.0.99) then
                         has_sigma_d = .true.
                         exit lp
                      endif
                   end do
                   !! - -
                   ax = - ax_list(:,i) - ax_list(:,j)
                   ax = ax/norm2(ax)
                   ! write(*,*) ax
                   !! dot with s0 should be zero -> s0 bisects axes i,j
                   do k = 1, nbas
                      if( op(k) .ne. 'S') cycle
                      if( n_int(k) .ne. 0) cycle
                      dotk = abs( dot_product(ax, ax_list(:,k)) )
                      call cross_prod( ax_list(:,k), cn_ax, cc)
                      cross=norm2(cc)
                      ! write(*,*) i,j,k, dotk
                      if( dotk .lt. 0.01 .and. cross.gt.0.99) then
                         has_sigma_d = .true.
                         exit lp
                      endif
                   end do

                endif
             end do
          endif
       end do lp

    end if



    !! check if has s2n
    do i = 1, nbas
       if( op(i) .ne. 'S' ) cycle
       if( n_int(i) .ne. 2*max_n_val) cycle
       has_s2n = .true.
       exit
    end do

    ! write(*,*) 'c2 perp:',c2_perp
    ! write(*,*) 'has sigma h:',has_sigma_h
    ! write(*,*) 'has sigma v:', has_sigma_v
    ! write(*,*) 'has s2n:',has_s2n
    ! write(*,*) 'has sigma d:', has_sigma_d

    !! PG is one of C# or D#
    if( has_cn )then
       !!
       if( c2_perp ) then
          !! is D#
          !! c2 are perp to cn
          write(pg,'(a1,g0)') 'D',max_n_val
          if( has_sigma_d) write(pg,'(a1,g0,a1)') 'D',max_n_val,'d'
          if( has_sigma_h) write(pg,'(a1,g0,a1)') 'D',max_n_val,'h'
       else
          !! is C#
          write(pg,'(a1,g0)') 'C',max_n_val
          if( has_s2n ) write(pg,'(a1,g0)') 'S',2*max_n_val
          if( has_sigma_v ) write(pg,'(a1,g0,a1)') 'C',max_n_val,'v'
          if( has_sigma_h ) write(pg,'(a1,g0,a1)') 'C',max_n_val,'h'
       end if
       !!
    else
       !!
       ! pg = 'C1'
       if( has_sigma ) pg = 'Cs'
       if( has_inversion ) pg = 'Ci'
    endif

100 continue

    if( verb) then
       write(*,*) " >> Report from get_pg:"
       write(*,*) 'has inversion:',has_inversion
       write(*,*) 'has sigma:', has_sigma
       write(*,*) 'has cn:', has_cn
       write(*,*) 'max_n_val:',max_n_val
       write(*,'(a14,3f9.4)') 'principal ax:',ax_list(:,max_n_loc)
       write(*,*) 'c2 perp:',c2_perp
       write(*,*) 'has sigma h:',has_sigma_h
       write(*,*) 'has sigma v:', has_sigma_v
       write(*,*) 'has s2n:',has_s2n
       write(*,*) 'has sigma d:', has_sigma_d
    end if



    !! some simple check on consistency of PG and number of Ops
    pg_err = check_pg_Nop( pg, nbas )
    ! write(*,*) pg
    ! if( pg_err .gt. 0 ) pg = '---'
    if( pg_err .gt. 0 ) then
       write(pg,'(a,a1)') trim(pg),'+'
    elseif( pg_err .lt. 0 )then
       write(pg,'(a,a1)') trim(pg),'-'
    endif

    ! write(*,*) 'PG is: ',pg

  end subroutine sofi_get_pg


  subroutine sofi_analmat( rmat, op, n, p, ax, angle )
    !!
    !! Analyse the input 3x3 matrix rmat, return the Schoenflies PG
    !! notation: "Op n^p", also give axis, angle of the operation.
    !! The angle is strictly positive value on output
    !! The +/- orientation of axis can be arbitrary, since it comes from eigenvector
    !!
    !! NOTE: two matrices M and M^T can produce the same result from this routine,
    !!   the relative orientation of operation can be lost (ambiguous).
    !!   To lift the ambiguity the whole list of operations needs to be processed,
    !!   see the routine sofi_unique_ax_angle for an attempt at this.
    !!
    use sofi_tools, only: diag, gcd_rec
    implicit none
    real, dimension(3,3), intent(in) :: rmat
    character(len=2),     intent(out) :: op   !! character for operation Id, I, C, S
    integer,              intent(out) :: n    !! the n for angle
    integer,              intent(out) :: p    !! power for C5^2 kinda stuff...
    real, dimension(3),   intent(out) :: ax   !! the axis
    real,                 intent(out) :: angle

    real, dimension(3,3) :: rdum
    real :: det, pi, search_eval, diff, diff_old, cosangl
    real, dimension(3) :: eigvals
    integer :: idx, i, j, nl, pl, gcd
    real :: mindiff

    op=''
    pi = 4.0*atan(1.0)

    !! working copy
    rdum(:,:) = rmat(:,:)

    !! determinant
    !! IRA routine
    call determinant3x3( rdum, det )

    !! diagonalise
    call diag( 3, rdum, eigvals, 1)

    ! write(*,*) 'diagonalised'
    ! do i = 1, 3
    !    write(*,'(3f9.4)') rdum(i,:)
    ! end do

    !! based on value of det, decide what to check
    if( det .gt. 0.5 ) then
       !! positive 1.0 determinant, matrix is rotation
       !! there is one eignevalue which is 1.0
       search_eval = 1.0
       ! write(*,*) "matrix is rotation"
       op(1:1)='C'
       !!
    elseif( det .lt. -0.5) then
       !! negative 1.0 determinant, matrix is (roto-)inversion.
       !! There is one eignevalue -1.0
       search_eval = -1.0
       ! write(*,*) "matrix is (roto-)inversion"
       op(1:1)='S'
    endif

    !! find requested eigenvalue (they are not ordered)
    diff_old = 99.9
    idx = 0
    do i = 1, 3
       diff = eigvals(i) - search_eval
       if( abs(diff) .lt. diff_old) then
          diff_old = abs(diff)
          idx = i
       end if
    end do

    !! ax is the corresponding eigenvector
    ax = rdum(:,idx)
    !! cosine of the angle is given by two other degenerate eigvals, select one.
    !! they are not sorted after diag...
    i = idx+1
    i = mod(i,3) + 1
    cosangl = eigvals(i)
    !! massage the eigenvalue a bit, since acos is defined strictly on [-1:1]
    if( cosangl .gt. 1.0 .or. cosangl .lt. -1.0) cosangl = real(nint(cosangl))
    angle = acos( cosangl ) / (2.0*pi)
    !! the resulting angle here is always positive

    ! write(*,*) cosangl, angle
    ! write(*,*) 'analmat angle:',angle*360

    !! Schoenflies notation: "OP n^p"
    !! primitive way to find n and p. There should be a better way for this...!
    !! loop over 12 => 1/12 seems to be smallest angle of rotation in a PG.
    !! Go beyond 12 ... ?
    n = 0
    p = 1
    nl = n
    pl = p
    mindiff=99.9
    do j = 1, 24
       diff = angle*j - nint(angle*j)
       ! write(*,*) j, diff
       if( abs(diff) .lt. mindiff ) then
          nl = j
          pl = nint(angle*j)
          ! exit
          mindiff = abs(diff)
       end if
    end do
    n = nl
    p = pl
    !! greatest common denominator
    gcd = gcd_rec( p, n)
    ! write(*,*) "gcd:",n, p, gcd
    n = n/gcd
    p = p/gcd


    if( angle .lt. 1e-3 ) n = 0
    if( p .eq. 0 ) p = 1
    !!
    !! C0 = Id :: rotation of angle 0 around any axis is identity
    if( op(1:1)=='C' .and. angle .lt. 1e-6 ) op='Id'
    !!
    !! S2 = I :: rotation 0.5 and reflection about any axis is inversion
    if( op(1:1)=='S' .and. abs(angle-0.5) .lt. 1e-6 ) op='I'

    ! write(*,'(a,x,i3,a3,2f9.4)') 'angle:',n, op, angle, (n-1.0/angle)

    !! detection of possible errors
    if( &
         op(1:1)=='C' .and. n==0 &
         ) then
       ! write(*,'(a,2x,a3,x,g0,x,g0)') '::: Error in sofi_analmat',op,n,p
       ! stop
    end if

  end subroutine sofi_analmat



  subroutine sofi_ext_Bfield( n_op, op_list, b_field )
    !! the effect of external B field on PG is to filter the list of Ops,
    !! only those which satisfy:
    !!
    !!    det( M ) M B = B
    !!
    !! are valid, where M is the symmOp matrix.
    !!
    !!==================================================
    !! see Eq(2) of:
    !! A. Pausch, M Gebele, W. Klopper, J. Chem. Phys. 155, 201101 (2021)
    !! https://doi.org/10.1063/5.0069859
    !!==================================================
    use sofi_tools, only: op_valid_ext_b
    implicit none
    integer,                           intent(inout) :: n_op
    real, dimension(1:3, 1:3, 1:n_op), intent(inout) :: op_list
    real, dimension(3),                intent(in) :: b_field

    integer :: i, n_op_new
    integer, dimension(n_op) :: which
    real, dimension(3,3) :: rmat
    real, dimension(3,3,n_op) :: op_new

    logical :: valid

    !! array to keep track of which ops survive the filter
    which(:) = 0

    n_op_new = 0
    op_new(:,:,:) = 0.0

    do i = 1, n_op
       rmat = op_list(:,:,i)
       !!
       !! validate rmat with B
       !!
       valid = op_valid_ext_b( rmat, b_field )
       !!
       if( valid ) then
       !! M satisfy the Eq.2
          which(i) = 1
          n_op_new = n_op_new + 1
          op_new(:,:,n_op_new) = rmat
       end if
    end do

    !! set output
    n_op = n_op_new
    op_list(:,:,:) = 0.0
    op_list(:,:,1:n_op) = op_new(:,:,1:n_op)

    if( n_op .lt. 1) then
       !! should not, at least Id should satisfy!
       write(*,*) 'error in sofi_ext_Bfield'
       stop
    endif

    return
  end subroutine sofi_ext_Bfield


  subroutine sofi_construct_operation( op, axis, angle, matrix )
    !! angle on input is in units 1/(2pi), e.g. angle=0.5 means half circle
    use sofi_tools, only: construct_rotation, construct_reflection
    implicit none
    character(len=2),     intent(in) :: op
    real, dimension(3),   intent(in) :: axis
    real,                 intent(in) :: angle
    real, dimension(3,3), intent(out) :: matrix

    real, dimension(3,3) :: tmp
    real, dimension(3) :: ax
    real :: an
    real, parameter :: twopi=8.0*atan(1.0)

    !! put angle into 2pi units
    an = angle*twopi
    !! normalize axis
    ax = axis/norm2(axis)
    !! initialize
    matrix(:,:) = 0.0
    !!
    select case( trim(adjustl(op)) )
    case( 'Id', 'E', 'id', 'e' )
       !! identity
       matrix(1,1) = 1.0
       matrix(2,2) = 1.0
       matrix(3,3) = 1.0
    case( 'I', 'i' )
       !! inversion
       matrix(1,1) = -1.0
       matrix(2,2) = -1.0
       matrix(3,3) = -1.0
    case( 'C', 'c' )
       !! rotation
       call construct_rotation( ax, an, matrix )
    case( 'S', 's' )
       !! rotoreflection
       call construct_rotation( ax, an, tmp )
       call construct_reflection( ax, matrix )
       matrix = matmul( matrix, tmp )
    case default
       write(*,*) "ERROR in sofi_construct_operation: unknown 'op':",trim(adjustl(op))
       stop
    end select

  end subroutine sofi_construct_operation

  subroutine sofi_unique_ax_angle( n_mat, mat_list, op_out, ax_out, angle_out )
    !! This routine is an attempt, should be taken with caution.
    !! Generate list of op, axis, angle which resolves ambiguity for
    !! conjugate operations (transpose)
    !!=======
    !! method:
    !!   analyse each matrix in default way, then reconstruct prescribed matrix from op, ax, angle.
    !!   Compare the result matrix with original, if equal, then angle is good, otherwise angle is negative.
    !!   NOTE: before doing this, the equivalent axes should be flipped to point in the same direction,
    !!         and this also changes the sign of angle. But it is not always that M and M^T have opposite
    !!         direction of axes, the +/- direction after diagonalization is random.

    use sofi_tools, only: matrix_distance
    implicit none
    integer,                    intent(in) :: n_mat
    real, dimension(3,3,n_mat), intent(in) :: mat_list
    character(len=2), dimension(n_mat), intent(out) :: op_out
    real, dimension(3,n_mat),           intent(out) :: ax_out
    real, dimension(n_mat),             intent(out) :: angle_out

    integer :: n, i, p, j
    real :: angle, dotp, angle_diff, dist, dist_neg
    character(len=2), dimension(n_mat) :: op_list
    real, dimension(3) :: ax
    real, dimension(3,3) :: rmat

    !! get all axes
    do i = 1, n_mat
       rmat = mat_list(:,:,i)
       call sofi_analmat( rmat, op_list(i), n, p, ax, angle )
       op_out(i) = op_list(i)
       ax_out(:,i) = ax
       angle_out(i) = angle
    end do


    !! flip equivalent axes, if ax flipped, flip also angle
    do i = 1, n_mat
       do j = 1, n_mat
          dotp = dot_product( ax_out(:,i), ax_out(:,j) )
          if( abs(dotp) .gt. 0.99 ) then
             !! axes are equivalent, flip
             ax_out(:,j) = dotp * ax_out(:,j)
             !! flip also angle
             angle_out(j) = dotp*angle_out(j)
          endif
       end do
       ! write(*,'(i4,x,3f9.4,4x,f9.4)') i, ax_out(:,i), angle_out(i)
    end do

    !! find operations that are still ambiguous
    do i = 1, n_mat
       ! write(*,'(a,i0,a4,3f9.4,4x,f9.4)') '>>',i, op_list(i), ax_out(:,i), angle_out(i)
       do j = i+1, n_mat
          !!
          dotp=dot_product( ax_out(:,i), ax_out(:,j))
          angle_diff = abs( angle_out(i) - angle_out(j) )
          !!
          ! write(*,*) i,j, dotp, angle_diff
          if( op_list(i) .eq. op_list(j) .and. &
              dotp .gt. 0.99 .and. &
              angle_diff .lt. 1e-3 ) then
             !!
             !! found ambiguous operations (i,j)
             !! could just randomly decide and flip either, but let's
             !! do it proper.
             !!
             ! write(*,*) 'ambiguous angle with j=',j
             ! write(*,*) 'imat'
             ! write(*,'(3f12.6)') mat_list(1,:,i)
             ! write(*,'(3f12.6)') mat_list(2,:,i)
             ! write(*,'(3f12.6)') mat_list(3,:,i)
             ! write(*,*) 'jmat'
             ! write(*,'(3f12.6)') mat_list(1,:,j)
             ! write(*,'(3f12.6)') mat_list(2,:,j)
             ! write(*,'(3f12.6)') mat_list(3,:,j)
             !!
             !! reconstruct rmat with this ax, +/-angle
             !! imat
             call sofi_construct_operation( op_list(i), ax_out(:,i), angle_out(i), rmat )
             call matrix_distance( mat_list(:,:,i), rmat, dist )
             call sofi_construct_operation( op_list(i), ax_out(:,i), -angle_out(i), rmat )
             call matrix_distance( mat_list(:,:,i), rmat, dist_neg )
             ! write(*,*) 'imat',dist, dist_neg
             if( dist_neg .lt. dist ) then
                !! angle i should be flipped
                ! write(*,*) 'flipping angle i'
                angle_out(i) = -angle_out(i)
             endif


             !! jmat
             call sofi_construct_operation( op_list(j), ax_out(:,j), angle_out(j), rmat )
             call matrix_distance( mat_list(:,:,j), rmat, dist )
             call sofi_construct_operation( op_list(j), ax_out(:,j), -angle_out(j), rmat )
             call matrix_distance( mat_list(:,:,j), rmat, dist_neg )
             ! write(*,*) 'jmat',dist, dist_neg
             if( dist_neg .lt. dist ) then
                !! angle j should be negative
                ! write(*,*) 'flipping angle j'
                angle_out(j) = -angle_out(j)
             endif

             !! the angles should now be different!
             if( abs(angle_out(i) - angle_out(j)) .lt. 0.01 ) then
                !! error!
                write(*,*) ">>>>! ERROR, angles still ambiguous!"
             endif


          endif
          !!
       end do
    end do

    !! put ax of Id or I ops to (1, 0, 0 )
    do i = 1, n_mat
       if( op_out(i) == "Id" ) ax_out(:,i) = (/ 1.0, 0.0, 0.0 /)
       if( op_out(i) == "I" ) ax_out(:,i) = (/ 1.0, 0.0, 0.0 /)
    end do



    !! output final list
    ! write(*,*) '----- final'
    ! do i = 1, n_mat
    !    write(*,'(i4,x,a,x,3f9.4,4x,f9.4)') i, op_out(i), ax_out(:,i), angle_out(i)
    ! end do



    !! reconstruct each from op, angle, ax
    !! reconstruct + angle

    !! reconstruct - angle
  end subroutine sofi_unique_ax_angle


  subroutine sofi_mat_combos( n_in, mat_in, n_out, mat_out )
    !! generate combos of matrices without a structure
    use sofi_tools, only: nmax, matrix_distance
    implicit none
    integer, intent(in) :: n_in
    real, dimension(3,3,n_in), intent(in) :: mat_in
    integer, intent(out) :: n_out
    real, dimension(3,3,nmax), intent(out) :: mat_out

    integer :: i, m, ii, j, nn, k
    real, dimension(3,3) :: rmat
    logical :: is_new
    real :: dd

    !! copy input matrices
    mat_out(:,:,1:n_in) = mat_in(:,:,:)

    nn = n_in
    mm: do m = 1, 10
       ii = nn
       do i = 1, nn
          do j = 1, nn
             !! construct combo
             rmat = matmul( mat_out(:,:,i), mat_out(:,:,j) )
             !! check if we already have it
             is_new = .true.
             pp: do k = 1, ii
                call matrix_distance( rmat, mat_out(:,:,k), dd )

                if( dd .lt. 0.73 ) then
                   is_new = .false.
                   exit pp
                endif

             end do pp
             !!
             if( is_new ) then
                !! add it
                ii = ii + 1
                mat_out(:,:,ii) = rmat
             endif
             !!
          end do
       end do
       !! if nothing new, exit
       if( nn .eq. ii) exit mm
       nn = ii
    end do mm

    n_out = nn

  end subroutine sofi_mat_combos