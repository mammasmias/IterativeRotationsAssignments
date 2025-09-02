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

!> @details
!!
!! Compute all data from SOFI in single routine.
!! This contains getting the list of symmetry operations,
!! and the associated permutations, Op symbols, n and p integers,
!! list of axis, angles, and dHausdorff. Also the PG name, and list of principal axes.
!!
!! The structure is assumed to be already shifted to desired origin on input.
!!
!! @warning
!!  This routine performs the combinations of found symmetry operations until
!!  group completeness. The ``sym_thr`` is disregarded for the combinations, therefore
!!  operations resulting as combinations of other operations can have a ``dHausdorff`` value higher
!!  than ``sym_thr``!
!!
!! @note
!!  All arrays have size ``nmax``, but the actual output values go only up to index ``nmat`` (or ``n_prin_ax``),
!!  values beyond that index can be random.
!!
!! @param[in] nat                        :: number of atoms
!! @param[in] typ(nat)      :: integer atomic types
!! @param[in] coords(3,nat) :: positions of atoms
!! @param[in] sym_thr       :: threshold for finding symmetries (not taken into account when making combinations)
!! @param[in] prescreen_ih  :: flag to check early-termination for Ih groups
!! @param[out] nmat         :: number of symmetries found
!! @param[out] mat_list(3,3,nmat)   :: symmetry matrices
!! @param[out] perm_list(nat, nmat) :: permutation of atoms after applying each symmetry
!! @param[out] op_list(nmat)        :: Character "Op" from the Schoenflies notation: Op n^p
!!                                 (E = identity, I = inversion, C = rotation, S = (roto-)reflection )
!! @param[out] n_list(nmat)       :: Schoenflies n value
!! @param[out] p_list(nmat)       :: Schoenflies p value
!! @param[out] ax_list(3,nmat)    :: axis of operation of each symmetry
!! @param[out] angle_list(nmat)   :: angle of each symmetry, in units of 1/2pi,
!!                               i.e. angle=0.333 is 1/3 of full circle, or 120 degrees
!! @param[out] dHausdorff_list(nmat)  :: max difference of atomic positions of before/after symm transformation
!! @param[out] pg           :: name of Point group, e.g. D6h
!! @param[out] n_prin_ax    :: number of equivalent principal axes
!! @param[out] prin_ax(3,n_prin_ax) :: list of equivalent principal axes.
!! @param[out] ierr         :: error value, zero on normal execution, negative otherwise
!! @returns nmat, mat_list, perm_list, op_list, n_list, p_list, ax_list, angle_list, dHausdorff_list, pg, prin_ax, ierr
!!
subroutine sofi_compute_all( nat, typ, coords, sym_thr, prescreen_ih, &
     nmat, mat_list, perm_list, &
     op_list, n_list, p_list, &
     ax_list, angle_list, dHausdorff_list, pg, n_prin_ax, prin_ax, &
     ierr )
  use ira_precision
  use sofi_tools, only: nmax
  use err_module
  implicit none
  !! ===== input
  integer(ip),                 intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ
  real(rp), dimension(3,nat),  intent(in) :: coords
  real(rp),                    intent(in) :: sym_thr
  logical,                 intent(in) :: prescreen_ih
  !! ===== output
  integer(ip),                       intent(out) :: nmat
  real(rp), dimension(3,3,nmax),     intent(out) :: mat_list
  integer(ip), dimension(nat, nmax), intent(out) :: perm_list
  character(len=1), dimension(nmax), intent(out) :: op_list
  integer(ip), dimension(nmax),      intent(out) :: n_list
  integer(ip), dimension(nmax),      intent(out) :: p_list
  real(rp), dimension(3, nmax),      intent(out) :: ax_list
  real(rp), dimension(nmax),         intent(out) :: angle_list
  real(rp), dimension(nmax),         intent(out) :: dHausdorff_list
  character(len=10),             intent(out) :: pg
  integer(ip),                       intent(out) :: n_prin_ax
  real(rp), dimension(3,nmax),       intent(out) :: prin_ax
  integer(ip),                       intent(out) :: ierr

  real(rp) :: dum
  integer(ip) :: i
  real(rp), dimension(3) :: rdum
  real(rp), dimension(3,3) :: rmat
  logical :: verb
  character(:), allocatable :: msg

  real(rp) :: dt


  !! get symmetries with sym_thr
  call sofi_get_symmops( nat, typ, coords, sym_thr, prescreen_ih, nmat, mat_list, ierr )
  if( ierr /= 0 ) then
     write(*,*) "at: ",__FILE__," line:",__LINE__
     return
  end if


  !! get combos, can produce symm above sym_thr
  call sofi_get_combos( nat, typ, coords, nmat, mat_list, ierr )
  if( ierr /= 0 ) then
     write(*,*) "at: ",__FILE__," line:",__LINE__
     return
  end if

  !! get ops, and unique angles and axes
  !! (potentially not needed due to convention on direction of axes)
  call sofi_unique_ax_angle( nmat, mat_list, op_list, ax_list, angle_list, ierr )


  if( ierr /= 0 ) then
     write(*,*) "at: ",__FILE__," line:",__LINE__
     return
  end if

  !! get permuations and dHausdorff
  call sofi_get_perm( nat, typ, coords, nmat, mat_list, perm_list, dHausdorff_list )


  !! get name of pg
  verb = .false.
  ! verb = .true.
  call sofi_get_pg( nmat, mat_list, pg, n_prin_ax, prin_ax, verb, ierr )
  if( ierr /= 0 ) then
     write(*,*) "at: ",__FILE__," line:",__LINE__
     return
  end if


  !! get the p, n values
  do i = 1, nmat
     rmat = mat_list(:,:,i)
     call sofi_analmat( rmat, op_list(i), n_list(i), p_list(i), rdum, dum, ierr )
     if( ierr /= 0 ) then
        write(*,*) "at: ",__FILE__," line:",__LINE__
        return
     end if

  end do

end subroutine sofi_compute_all


!> @details
!! Wrapper routine for directly getting only the PG name from a structure, nothing else.
!! The structure is shifted to geometrical center in this routine!
!!
!! @warning
!!  This routine performs the combinations of found symmetry operations until
!!  group completeness. The ``sym_thr`` is disregarded for the combinations, therefore
!!  operations resulting as combinations of other operations can have a ``dHausdorff`` value higher
!!  than ``sym_thr``!
!!
!! @param[in] nat :: number of atoms
!! @param[in] typ_in(nat) :: atomic species
!! @param[in] coords_in(3,nat) :: atomic positions
!! @param[in] sym_thr :: symmetry threshold
!! @param[out] pg :: point group
!! @param[in] verb :: flag for verbose output
!!
subroutine sofi_struc_pg( nat, typ_in, coords_in, sym_thr, pg, verb )
  use ira_precision
  use sofi_tools, only: nmax
  implicit none
  integer(ip),                 intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ_in
  real(rp), dimension(3,nat),  intent(in) :: coords_in
  real(rp),                    intent(in) :: sym_thr
  character(len=10),       intent(out) :: pg
  logical,                 intent(in) :: verb
  ! local
  integer(ip) :: n_op
  integer(ip), allocatable :: perm_list(:,:)
  real(rp), allocatable :: op_list(:,:,:)
  real(rp), dimension(3) :: gc
  real(rp), allocatable :: prin_ax(:,:)
  integer(ip) :: n_prin_ax
  integer(ip) :: i, ierr
  integer(ip), dimension(nat) :: typ
  real(rp), dimension(3,nat) :: coords
  logical :: prescreen_ih

  pg = ""

  allocate( op_list(1:3,1:3,1:nmax) )
  allocate( perm_list(1:nat, 1:nmax) )
  allocate( prin_ax(1:3,1:nmax) )

  !! working copies
  typ = typ_in
  coords = coords_in

  !! recenter to gc
  gc = sum( coords(:,:),2)/nat
  do i = 1, nat
     coords(:,i) = coords(:,i) - gc
  end do

  !! get list of symm operations
  prescreen_ih = .false.
  ! prescreen_ih = .true.
  call sofi_get_symmops( nat, typ, coords, sym_thr, prescreen_ih, n_op, op_list, ierr )
  ! call sofi_get_symmops( nat, typ, coords, sym_thr, n_op, op_list, perm_list )
  if( ierr /= 0 ) return

  !! get PG name
  ! call sofi_get_pg( n_op, op_list, pg, verb )

  !! combos with higher sym_thr: find the missing Ops in case of '+' or '-' PGs
  call sofi_get_combos( nat, typ, coords, n_op, op_list, ierr )
  if( ierr /= 0 ) return
  ! call sofi_get_combos( nat, typ, coords, n_op, op_list, perm_list )

  !! get new PG name
  call sofi_get_pg( n_op, op_list, pg, n_prin_ax, prin_ax, verb, ierr )
  if( ierr /= 0 ) return

  deallocate( op_list, perm_list, prin_ax )

end subroutine sofi_struc_pg



!> @details
!!
!! The main SOFI routine
!!======================
!!
!! Find the list of symmetry operations of the atomic structure, with a threshold ``sym_thr``.
!! The symmetry operations are in form of 3x3 matrices.
!!
!! When distortions are present in the atomic structure, the distance between the original
!! structure \f$A\f$ and the symmetry-transformed structure \f$\theta A\f$ is some small
!! nonzero value. This distance is computed as largest distance of any atom after the transformation
!! from the position of the corresponding atom in the original structure (Hausdorff distance).
!! The threshold ``sym_thr`` gives the upper limit for this distance value.
!!
!! @note
!!  The output ``op_list`` is assumed to be allocated by the caller, and has the size [3,3,nmax],
!!  where ``nmax=200`` by default, but can be adjusted in sofi_tools.f90
!!
!! The operations from ``op_list`` can be applied by matmul(), for example the N-th entry:
!!
!!~~~~~~~~~~~~~~{.f90}
!!    theta(3,3) = op_list(:, :, N)
!!    do i = 1, natoms
!!       coords_transformed(:,i) = matmul( theta, coords_original(:,i) )
!!    end do
!!~~~~~~~~~~~~~~
!!
!! @param[in] nat              :: number of atoms
!! @param[in] typ_in(nat)      :: integer atomic types
!! @param[in] coords_in(3,nat) :: positions of atoms
!! @param[in] sym_thr        :: threshold for finding symmetries, in terms of Hausdorff distance
!! @param[in] prescreen_ih   :: flag to check early termniation for Ih
!! @param[out] n_so          :: number of found symmetry operations
!! @param[out] op_list       :: the list of operations
!! @param[out] ierr          :: error value, negative on error, zero otherwise
!! @returns n_so, op_list, ierr
!!
subroutine sofi_get_symmops( nat, typ_in, coords_in, sym_thr, prescreen_ih, n_so, op_list, ierr )
  use ira_precision
  use sofi_tools, only: nmax, m_thr, construct_reflection
  use err_module
#ifdef DEBUG
  use timer
#endif
  implicit none
  integer(ip),                 intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ_in
  real(rp), dimension(3,nat),  intent(in) :: coords_in
  real(rp),                    intent(in) :: sym_thr
  logical,                 intent(in) :: prescreen_ih
  !!
  integer(ip),                           intent(out) :: n_so
  real(rp), dimension(3,3,nmax),         intent(out) :: op_list
  integer(ip),                           intent(out) :: ierr

  !! local
  integer(ip), dimension(nat) :: typ
  real(rp), dimension(3,nat) :: coords
  real(rp), allocatable :: d_o(:,:)
  integer(ip) :: i, j, k, l, m, mm
  real(rp), dimension(3,3) :: theta, beta, gamma
  logical :: fail1, fail_beta, is_collinear
  real(rp) :: dd, d_i, d_j, dh
  integer(ip) :: nbas, nn
  integer(ip) :: ti, tj
  real(rp) :: small_norm
  !! state
  logical :: has_sigma, has_inversion, terminate, success, prescreen
  integer(ip) :: n_old, n_u_c3
  real(rp) :: u_c3(3,5)
  real(rp) :: ax(3)
  real(rp) :: rmin, rmax
#ifdef DEBUG
  type( local_timer ) :: tm
#endif

#ifdef DEBUG
  call tm%tag(1, "symmops")
  call tm%start(1)
#endif


  ierr = 0
  prescreen = prescreen_ih

  if( size(op_list,3) .ne. nmax ) then
     ierr = ERR_SIZE_NMAX
     write(*,*) 'the op_list on input should have maxsize:',nmax
     write(*,*) "origin at: ",__FILE__," line:",__LINE__
     return
  end if

  !! edge cases for nat=2 or nat=1
  if( nat .lt. 3 ) then
     !! add identity

     !! if nat==2: check mirror

     !! return
  end if

  !! zero the output data
  n_so = 0
  ! op_list(:,:,:) = 0.0_rp

  !! set local copy
  typ(:) = typ_in(:)
  coords(:,:) = coords_in(:,:)

  !! sort by d
  rmin = 999.9_rp
  rmax = 0.0_rp
  allocate( d_o(1:2,1:nat))
  do i = 1, nat
     d_o(1,i) = norm2(coords(:,i))
     d_o(2,i) = real(i, rp)
     rmin = min(rmin, d_o(1,i))
     rmax = max(rmax, d_o(1,i))
  end do
  !! heuristic "detection" of fullerene ... impose prescreen
  if( rmax - rmin < 0.5_rp ) prescreen = .true.

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
  dd = 99.0_rp
  do i = 1, nat
     do j = i+1, nat
        dd = min( dd, norm2(coords(:,i) - coords(:,j)))
     end do
  end do

  !! initialize the state
  n_u_c3 = 0
  u_c3(:,:) = 0.0_rp
  has_sigma = .false.
  has_inversion = .false.

  !! add id
  theta(:,:) = 0.0_rp
  do i = 1, 3
     theta(i,i) = 1.0_rp
  end do

  call try_sofi( theta, nat, typ, coords, sym_thr, dd, n_so, op_list, dh, 0.5_rp, success, ierr )
  if( ierr /= 0 ) then
     write(*,*) "at: ",__FILE__," line:",__LINE__
     return
  end if


  !! try inversion
  theta(:,:) = 0.0_rp
  do i = 1, 3
     theta(i,i) = -1.0_rp
  end do

  call try_sofi( theta, nat, typ, coords, sym_thr, dd, n_so, op_list, dh, 0.5_rp, success, ierr )
  if( ierr /= 0 ) then
     write(*,*) "at: ",__FILE__," line:",__LINE__
     return
  end if


  !! edge case: check if structure is collinear
  ! is_collinear = check_collinear( nat, coords, ax )
  call sofi_check_collinear( nat, coords, is_collinear, ax )
  if( is_collinear ) then
     ! write(*,*) "SOFI:: structure is collinear"
     !! construct reflection on the ax
     call construct_reflection( ax, theta )
     !! try... should be detected as I above, but test again anyway.
     call try_sofi( theta, nat, typ, coords, sym_thr, dd, n_so, op_list, dh, 0.5_rp, success, ierr )
     if( ierr /= 0 ) then
        write(*,*) "at: ",__FILE__," line:",__LINE__
        return
     end if
     return
  end if


  !! set beta
  small_norm = max(sym_thr, 1e-1_rp)

  !! fail1 happens when vectors are collinear, i.e. orthonormal basis could not be made
  fail1 = .true.
  !! fail_beta is .true. when no beta is generated
  fail_beta = .true.
  do i = 1, nat
     if( norm2( coords(:,i)) .lt. small_norm) cycle
     do j = 1, nat
        if( norm2( coords(:,j)) .lt. small_norm) cycle
        !!
        !! call to IRA lib
        call set_orthonorm_bas( coords(:,i), coords(:,j), beta, fail1 )
        if( fail1 ) cycle
        !!
        fail_beta = .false.
        exit
     end do
     exit
  end do

  ! write(*,*) "fail_beta", fail_beta
  ! write(*,*) i,j
  if( fail1 .or. i .gt. nat .or. j .gt. nat ) then
     write(*,*) repeat('%',40)
     write(*,*) "ERROR: cannot set beta. Structure not properly shifted?"
     write(*,*) repeat('%',40)
     write(*,*) "small norm:",small_norm
     write(*,*) nat
     write(*,*) "i=",i,"j=",j
     do i = 1, nat
        write(*,*) typ(i), coords(:,i)
     end do
     ierr = ERR_BETA
     write(*,*) "origin at: ",__FILE__," line:",__LINE__
     return
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
  kloop_: do k = 1, nat
     if( abs( d_i - norm2(coords(:,k))) .gt. 1.1_rp*sym_thr ) cycle
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
#ifdef DEBUG
           call tm%tag(2, "try_sofi")
           call tm%start(2)
#endif

           call try_sofi( theta, nat, typ, coords, sym_thr, dd, n_so, op_list, dh, m_thr, success, ierr )

#ifdef DEBUG
           call tm%stop(2)
#endif

           if( ierr /= 0 ) then
              write(*,*) "at: ",__FILE__," line:",__LINE__
              return
           end if

           !! mat is not successful, go to next
           ! if( .not. success ) cycle

           if( success .and. prescreen ) then
              n_old = n_so

#ifdef DEBUG
              call tm%tag(3, "get_combos")
              call tm%start(3)
#endif

              call sofi_get_combos( nat, typ, coords, n_so, op_list, ierr )
              if( ierr /= 0 ) then
                 write(*,*) "at: ",__FILE__," line:",__LINE__
                 return
              end if


#ifdef DEBUG
              call tm%stop(3)
#endif


              call state_update( n_old, n_so, op_list, has_sigma, has_inversion, n_u_c3, u_c3, terminate )
              ! write(*,*) "current state print:"
              ! write(*,*) "has_sigma", has_sigma
              ! write(*,*) "has_inversion",has_inversion
              ! write(*,*) "n_u_c3", n_u_c3
              ! write(*,*) "terminate", terminate
              ! do i = 1, n_u_c3
              !    write(*,"(3f9.4)") u_c3(:,i)
              ! end do
              if( terminate ) exit kloop_
           end if

           ! if( success ) then
           !    !! TODO
           !    !! if op = OP_ROT, and angle is on edge of m_thr, generate the half-angle and try
           ! end if


           ! write(*,*) 'op exiting try_sofi:'
           ! call sofi_analmat( theta, ... )
           ! write(*,'(2i4,3f9.5,x,f9.4,a5,f9.4)') k,l,ax1,d1,'dh:',dh
        end do
     end do
  end do kloop_
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

#ifdef DEBUG
  call tm%stop(1)
  write(*,*) "symmops timer;"
  call tm%print()
#endif


  ! write(*,*) 'exiting get_symmops'
  deallocate( d_o )

end subroutine sofi_get_symmops



!> @details
!! Obtain the permutations of atoms associated to each of the symmetry operations in
!! the list. Also compute the dHausdorff of each operation.
!!
!! This consists of computing P by cshda, permuting, applying SVD, and finally computing dHausdorff.
!!
!! Applying the permutation associated to N-th symmetry operation:
!!~~~~~~~~~~~~~~{.f90}
!!    coords(:, perm_list(:, N) ) = coords(:,:)
!!~~~~~~~~~~~~~~
!!
!! @param[in] nat :: number of atoms
!! @param[in] typ(nat) :: atomic types
!! @param[in] coords(3,nat) :: atomic positions
!! @param[in] nbas :: number of operations in the list
!! @param[in] bas_list(3,3,nbas) :: list of symmetry operations
!! @param[out] perm_list(nat,nbas) :: list of permutations
!! @param[out] dHausdorff_list(nbas) :: list of dHausdorff
!! @returns perm_list, dHausdorff_list
!!
subroutine sofi_get_perm( nat, typ, coords, nbas, bas_list, perm_list, dHausdorff_list )
  use ira_precision
  implicit none
  integer(ip),                   intent(in) :: nat
  integer(ip), dimension(nat),   intent(in) :: typ
  real(rp), dimension(3,nat),    intent(in) :: coords
  integer(ip),                   intent(in) :: nbas
  real(rp), dimension(3,3,nbas), intent(in) :: bas_list
  integer(ip), dimension(nat, nbas), intent(out) :: perm_list
  real(rp), dimension(nbas),     intent(out) :: dHausdorff_list

  integer(ip) :: i, j, ierr
  real(rp), dimension(3,3) :: rmat
  integer(ip), dimension(nat) :: t_local
  real(rp), dimension(3,nat) :: c_local
  integer(ip), dimension(nat) :: found
  real(rp), dimension(nat) :: dists

  real(rp) :: svd_rot(3,3), svd_tr(3), rdum(3)

  do i = 1, nbas
     !! get matrix from list
     rmat = bas_list(:,:,i)
     !! transform struc
     do j = 1, nat
        t_local(j) = typ(j)
        c_local(:,j) = matmul(rmat, coords(:,j) )
     end do
     !! compute cshda with high thr (accept anything)
     call cshda( nat, typ, coords, nat, t_local, c_local, 99.9_rp, found, dists )
     ! if( any(found) .eq. 0 ) then
     !    write(*,*) "!! ERROR IN: sofi_get_perm, failed cshda?"
     !    stop
     ! endif

     !! correct with svd tr, rot should be diag(1,1,1), tr could be something small
     call svd_forcerot( nat, typ, coords, nat, t_local(found), c_local(:,found), svd_rot, svd_tr, ierr )
     if( ierr /= 0 ) return
     ! write(*,*) "rot after svd"
     ! write(*,*) svd_rot(1,:)
     ! write(*,*) svd_rot(2,:)
     ! write(*,*) svd_rot(3,:)
     ! write(*,*) "tr after svd"
     ! write(*,*) svd_tr

     !! apply svd
     do j = 1, nat
        c_local(:,j) = matmul( svd_rot, c_local(:,j) ) + svd_tr
     end do

     !! recompute dists, now should be good
     dists(:) = 0.0_rp
     do j = 1, nat
        rdum = coords(:,j) - c_local(:,found(j))
        dists(j) = norm2( rdum )
     end do



     perm_list(:,i) = found(1:nat)
     dHausdorff_list(i) = maxval( dists(1:nat) )
     ! write(*,*) i, maxval(dists(1:nat))
     ! write(*,*) rmat(1,:)
     ! write(*,*) rmat(2,:)
     ! write(*,*) rmat(3,:)

  end do

end subroutine sofi_get_perm



!> @details
!! find all combos with high sym_thr: this will accept anything.
!! However, some strange matrices cannot be constructed from combos of
!! matrices which are already found as symmetries with sym_thr value.
!!
subroutine sofi_get_combos( nat, typ, coords, nbas, bas_list, ierr )
  use ira_precision
  use sofi_tools, only: m_thr, nmax
  implicit none
  integer(ip), intent(in) :: nat
  integer(ip), dimension(nat), intent(in) :: typ
  real(rp), dimension(3, nat), intent(in) :: coords
  integer(ip), intent(inout) :: nbas
  real(rp), dimension(3, 3, nmax), intent(inout) :: bas_list
  integer(ip), intent(out) :: ierr
  ! integer(ip), dimension(nat, nmax), intent(inout) :: perm_list

  integer(ip) :: m, i, j, ii
  real(rp), dimension(3,3) :: theta
  real(rp) :: dh, dd, sym_thr
  logical :: success

  !! hard-coded thresholds (accept anything)
  dd = 10.0_rp
  sym_thr = 5.0_rp

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
           call try_sofi( theta, nat, typ, coords, sym_thr, dd, ii, bas_list, dh, m_thr, success, ierr )
           if( ierr /= 0 ) then
              write(*,*) "at: ",__FILE__," line:",__LINE__
              return
           end if
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


!> @details
!!
!! Try if the matrix ``theta`` gives dH within 5*sym_thr, then send it to refine.
!! If matrix is valid and new, add it to op_list
!!
!! ``nbas`` and ``op_list`` on input contain the currently known symmetry operations,
!! and on output the same info updated accordingly if ``theta`` is new or not.
!!
!! @param[inout] theta :: 3x3 trial matrix
!! @param[in] nat :: number of atoms
!! @param[in] typ_in(nat) :: atomic species
!! @param[in] coords_in(3,nat) :: atomic positions
!! @param[in] sym_thr :: symmery threshold in terms of dH
!! @param[in] dd :: thr for first cshda, "smallest atom-atom distance"; also used for detemining "good-enough" trial matrix theta
!! @param[inout] nbas :: number of operations
!! @param[inout] op_list(3,3,nbas) :: list of operations
!! @param[out] dh :: Hausdorff distance value
!! @param[in] m_thr :: threshold for matrix_distance checking if matrix is new or not
!! @param[out] ierr :: error value, negative on error, zero otherwise
!! @returns theta, nbas, op_list, dh
!!
subroutine try_sofi( theta, nat, typ_in, coords_in, sym_thr, dd, nbas, op_list, dh, m_thr, success, ierr )
  use ira_precision
  use sofi_tools, only: nmax, epsilon
  implicit none
  real(rp), dimension(3,3),          intent(inout) :: theta
  integer(ip),                       intent(in) :: nat
  integer(ip), dimension(nat),       intent(in) :: typ_in
  real(rp), dimension(3,nat),        intent(in) :: coords_in
  real(rp),                          intent(in) :: sym_thr
  real(rp),                          intent(in) :: dd
  integer(ip),                       intent(inout) :: nbas
  real(rp), dimension(3,3,nmax),     intent(inout) :: op_list
  real(rp),                          intent(out) :: dh
  real(rp),                          intent(in) :: m_thr
  logical,                       intent(out) :: success
  integer(ip),                       intent(out) :: ierr

  integer(ip), dimension(nat) :: typ
  real(rp), dimension(3,nat) :: coords
  integer(ip) :: i
  real(rp), dimension(nat) :: dists
  integer(ip), dimension(nat) :: found, perm
  logical :: not_crazy, is_new
  logical :: is_valid, do_refine
  real(rp) :: act_thr

  ! write(*,*) "try_sofi received theta:"
  ! write(*,'(3f9.4)') theta

  ierr = 0
  success = .false.

  !! local copy
  typ = typ_in
  coords = coords_in

  dh = 99.9_rp
  is_valid = .false.

  !! apply theta
  do i = 1, nat
     coords(:,i) = matmul( theta, coords(:,i) )
  end do

  !! actual thr used to decide if theta is "close enough" to refine.
  !! In principle anything beyond 0.5*dd allows for permuttaions of atoms,
  !! but due to distortions in positions, need to allow a bit more, since
  !! permutations could be "partially correct", and refine can correct that.
  act_thr = 0.6_rp*dd + epsilon

  ! write(*,*) nat*2
  ! write(*,*)
  ! do i = 1, nat
  !    write(*,*) 1, coords_in(:,i)
  !    write(*,*) 2, coords(:,i)
  ! end do


  !! check if theta is already known
  call is_new_sofi( theta, nbas, op_list, m_thr, is_new )


  ! write(*,*) 'is new',is_new
  if( .not. is_new ) then
     return
  end if


  !! call to IRA lib
  !! get first cshda
  call cshda( nat, typ_in, coords_in, &
       nat, typ, coords, &
       act_thr, found, dists )
  dh = maxval(dists,1)
  ! write(*,'(x,a,x,f9.4)') 'initial dh',dh

  ! write(*,*) dh, 5.0*sym_thr, 1.1*dd
  ! do i = 1, nat
  !    write(*,*) i, found(i), dists(i)
  ! end do
  ! write(*,*) '1st dh:',dh, sum(dists**2)
  ! write(*,*) "initial", dh, act_thr

  !! add another epsilon to act_thr, for numerics...
  if( dh .gt. act_thr+epsilon ) then
     ! write(*,*) "return:: dh is too big"
     return
  end if



  not_crazy = .true.
  !! cshda returned no assignments
  if( any(found .eq. 0)) not_crazy = .false.

  ! write(*,*) 'not_crazy',not_crazy
  if( .not. not_crazy ) then
     ! write(*,*) "return:: crazy csdha"
     return
  end if

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

  do_refine = .true.
  !! if first cshda already super good, no need to refine
  if( dh .lt. min(sym_thr*1e-3_rp, epsilon) ) then
     do_refine = .false.
     is_valid = .true.
  end if



  if( do_refine ) then

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
        call add_sofi( nat, theta, nbas, op_list, ierr )
        if( ierr /= 0 ) then
           write(*,*) "at: ",__FILE__," line:",__LINE__
           return
        end if
        success = .true.
     end if
  end if

  ! write(*,*) "out theta:"
  ! write(*,'(3f9.4)') theta
  ! write(*,*) '::: exit try_sofi', success, dh

end subroutine try_sofi


!> @details
!! Refine the ``theta`` matrix, which is on input "close-to" some real symmetry operation,
!! and on output should be the matrix of symmetry operation that minimizes SVD.
!!
!! Do few steps of ICP-like algo... until n steps, or until no new permutations are found.
!!
!! @note
!!  assume atoms are already permuted at input
!!
!! *_ref is static, original structure, *_in is the structure transformed by theta and P on input
!!
!! @param[in] nat_ref :: number of atoms
!! @param[in] typ_ref(nat_ref) :: atomic types
!! @param[in] coords_ref(3,nat_ref) :: atomic positions
!! @param[in] nat :: number of atoms
!! @param[in] typ_in(nat) :: atomic types
!! @param[in] coords_in(3,nat) :: atomic positions
!! @param[inout] theta(3,3) :: trial matrix
!! @param[inout] dh :: Hausdorff value
!! @param[out] perm(nat) :: final permutation
!! @returns theta, dh, perm
!!
subroutine refine_sofi( nat_ref, typ_ref, coords_ref, nat, typ_in, coords_in, theta, dh, perm )
  use ira_precision
  implicit none
  integer(ip),                     intent(in) :: nat_ref
  integer(ip), dimension(nat_ref), intent(in) :: typ_ref
  real(rp), dimension(3,nat_ref),  intent(in) :: coords_ref
  integer(ip),                     intent(in) :: nat
  integer(ip), dimension(nat),     intent(in) :: typ_in
  real(rp), dimension(3,nat),      intent(in) :: coords_in
  real(rp), dimension(3,3),        intent(inout) :: theta
  real(rp),                        intent(inout) :: dh
  integer(ip), dimension(nat),     intent(out) :: perm

  integer(ip) :: i, n, k
  integer(ip), dimension(nat) :: found, found_new
  real(rp), dimension(nat) :: dists
  integer(ip), dimension(nat) :: typ
  real(rp), dimension(3,nat) :: coords
  real(rp), dimension(3,3) :: svd_r
  real(rp), dimension(3) :: svd_t  !, ax

  integer(ip) :: ierr
  ! real(rp) :: dd, det1

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
          nat, typ, coords, svd_r, svd_t, ierr )
     if( ierr /= 0 ) then
        write(*,*) "SVD failure"
        write(*,*) "origin at: ",__FILE__," line:",__LINE__
        return
     end if


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
          99.0_rp, found_new, dists )

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

  ! write(*,*) "refined in",i

end subroutine refine_sofi


!> @details
!! Check if rmat is new in op_list.
!! The choice is made by computing matrix_distance between rmat and all matrices in op_list,
!! if the distance is below m_thr for any matrix, then rmat is considered already known.
!!
!! @param[in] rmat(3,3) :: trial matrix
!! @param[in] nbas :: number of operations in the list
!! @param[in] op_list(3,3,nmax) :: list of symmetry operations to check
!! @param[in] m_thr :: threshold for matrix_distance
!! @param[out] is_new :: true if rmat is new, and should be added to op_list
!! @returns is_new
!!
subroutine is_new_sofi( rmat, nbas, op_list, m_thr, is_new )
  use ira_precision
  use sofi_tools, only: nmax, matrix_distance
  implicit none
  real(rp), dimension(3,3),      intent(in) :: rmat
  integer(ip),                   intent(in) :: nbas
  real(rp), dimension(3,3,nmax), intent(in) :: op_list
  real(rp),                      intent(in) :: m_thr
  logical,                   intent(out) :: is_new

  integer(ip) :: i
  logical :: eq_rmat
  real(rp) :: dd
  real(rp) :: det1, det2, ddet
  real(rp) :: matrix_thr

  !! Currently it is decided to compare rmat matrix to all the previously found
  !! matrices to decide whether rmat is new or not.
  !! Alternatively, this decision could be to compare permutation 'perm' to all the
  !! permutations in 'perm_list', since each symmOp should have a unique permutation.
  !! Or not? Mirror triangle over the plane gives same permutation as orig...

  !! This could not be decided only based on det and tr of matrices,
  !! since M and M^T have identical values for det and tr.
  !! The axes and angles of each operation are not known at this point.

  !! threshold on equivalence of matrices, higher value can speed up overall calc,
  !! but too high can think different matrices are equal.
  !! Value around 0.5 corresponds to a
  !! rotation of matrix by cca 1/18*2pi, which seems generally ok.
  !! Should this be let as input?! Done.
  ! matrix_thr = 0.5
  matrix_thr = m_thr

  !! initial values
  is_new = .true.
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
        if( ddet .lt. 1e-3_rp ) then
           eq_rmat = .true.
        endif

        ! write(*,*) 'rmat'
        ! write(*,'(3f9.4)') rmat
        ! write(*,*) 'op_list(:,:,i)'
        ! write(*,'(3f9.4)') op_list(:,:,i)
     endif

     is_new = .not.eq_rmat

     if( .not. is_new ) return
  end do

  ! write(*,*) "is new:", is_new, dd, matrix_thr
end subroutine is_new_sofi


!> @details
!! Add rmat into mat_list, increase nbas by 1.
subroutine add_sofi( nat, rmat, nbas, mat_list, ierr )
  use ira_precision
  use sofi_tools, only: nmax
  use err_module
  implicit none
  integer(ip),                      intent(in) :: nat
  real(rp), dimension(3,3),         intent(in) :: rmat
  integer(ip),                      intent(inout) :: nbas
  real(rp), dimension(3,3,nmax),    intent(inout) :: mat_list
  integer(ip),                      intent(out) :: ierr


  !! increment nbas
  nbas = nbas + 1
  if( nbas .gt. nmax ) then
     ierr = ERR_LIST_TOO_SMALL
     write(*,*) "subroutine add_sofi::: ERROR ADDING RMAT, LIST TOO SMALL"
     write(*,*) "origin at: "__FILE__," line:",__LINE__
     return
  end if

  mat_list(:,:,nbas) = rmat
  ierr = 0

end subroutine add_sofi


!> @details
!! flowchart to determine PG notation from op_list
!! online: https://symotter.org/assets/flowchart.pdf
!!
!! @note
!!  The principal axis is found as axis of the proper rotation operation (C symbol),
!!  with the largest n value. In case of multiple equivalent principal axes,
!!  they are found as all unique axes with those properties.
!!
!! @note
!!  In cases when the list of operations is incomplete for a certain PG,
!!  this routine will append either + or - sign to the PG tag, signifying
!!  either (+) a group with additional operations, or (-) a PG with some
!!  missing operations. These are deduced from the expected total number of operations
!!  of a PG. Such case typically happens when atomic structure is highly disordered.
!!
!! @param[in] nbas :: number of symmetry operations
!! @param[in] op_list(3,3,nbas) :: the list of symmetry operations (3x3 matrices)
!! @param[out] pg(10) :: the point group tag
!! @param[out] n_prin_ax :: number of equivalent principal axes
!! @param[out] prin_ax(3,n_prin_ax) :: list of equivalent principal axes of PG
!! @param[in] verb :: flag for verbosity
!! @param[out] ierr :: error value, zero on normal execution, negative otherise
!! @returns pg, prin_ax
!!
subroutine sofi_get_pg( nbas, op_list, pg, n_prin_ax, prin_ax, verb, ierr )
  !!
  use ira_precision
  use sofi_tools
  implicit none
  integer(ip),                   intent(in) :: nbas
  real(rp), dimension(3,3,nbas), intent(in) :: op_list
  character(len=10),         intent(out) :: pg
  integer(ip),                   intent(out) :: n_prin_ax
  real(rp), dimension(3,nbas),   intent(out) :: prin_ax
  logical,                   intent(in) :: verb
  integer(ip),                   intent(out) :: ierr

  real(rp), dimension(3) :: ax, cn_ax
  real(rp), dimension(3,3) :: rmat
  real(rp) :: angle, dd, dot, cross
  integer(ip) :: i, j
  character(len=1), dimension(nbas) :: op
  integer(ip), dimension(nbas) :: n_int, power
  real(rp), dimension(3,nbas) :: ax_list
  integer(ip), dimension(nbas) :: skip_ax
  integer(ip), dimension(nbas) :: multip_ax !! total number of operations on the axis of this op
  integer(ip) :: max_n_val, max_n_loc
  integer(ip) :: nax, nr_n
  real(rp), dimension(5,nbas) :: uniq_ax  !! 4th dim is largest Cn on this ax, 5th is counter on op
  logical :: has_inversion, has_sigma, has_cn, has_sigma_h, has_sigma_v, has_s2n, c2_perp, has_sigma_d
  integer(ip) :: n_c2, n_in_plane, ndum
  real(rp), dimension(3) :: cc
  integer(ip) :: k, pg_err
  real(rp) :: dotj, dotk

  real(rp), dimension(nbas) :: angle_list
  integer(ip) :: cn_multip, cn_val
  logical :: isnew

  if( verb ) then
     write(*,*) repeat('=',20)
     write(*,*) "number of SymmOps entering get_pg:",nbas
  endif

  ierr = 0

  !! start from below
  pg = 'C1'
  if( nbas .eq. 1 ) then
     !! Structure is C1 (no symmetry except identity).
     !! No need to do anything in this routine.
     !!
     !! output single principal axis. Maybe rather set to zero?
     n_prin_ax = 1
     prin_ax(:,1) = (/1.0_rp, 0.0_rp, 0.0_rp/)
     !!
     return
  end if


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
     call sofi_analmat( rmat, op(i), n_int(i), power(i), ax_list(:,i), angle, ierr )
     if( ierr /= 0 ) then
        write(*,*) "at: ",__FILE__," line:",__LINE__
        return
     end if
     angle_list(i) = angle
     if(verb) write(*,'(i2,a3,1x,i0,"^",i0,2x,"axis:",1x,3(f9.6,1x),2x,"angle:",1x,g0.4)') &
          i, op(i), n_int(i), power(i), ax_list(:,i), angle
  end do

  uniq_ax(:,:) = 0.0_rp
  nax = 0
  multip_ax(:) = 0
  skip_ax(:) = 0


  !! some basics
  has_inversion = find_inversion( nbas, op )
  has_sigma = find_sigma( nbas, op, n_int )
  has_cn = find_cn( nbas, op, n_int )


  !! group the Ops of the list by axis
  do i = 1, nbas
     if( op(i) .eq. OP_INVERSION ) cycle
     if( op(i) .eq. OP_IDENTITY ) cycle
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
        !! skip E and I
        if( op(j) == OP_IDENTITY )cycle
        if( op(j) == OP_INVERSION )cycle
        !!
        if( abs(dot_product(ax, ax_list(:,j))) .gt. 0.999_rp ) then
           !! is same ax
           if(verb) write(*,'(a3,g0,a1,g0,1x,f7.4)') op(j), n_int(j),'^', power(j), angle_list(j)
           skip_ax(j) = 1
           !! keep maximal n value of C operations
           if( op(j) == OP_PROP_ROT) max_n_val = max( max_n_val, n_int(j))
           !! count how many ops have this ax
           uniq_ax(5,nax) = uniq_ax(5,nax) + 1.0_rp
        end if
     end do
     uniq_ax(4,nax) = real(max_n_val, rp)
  end do


  !! put multiplicity of ax to each op
  do i = 1, nbas
     !! find which ax from uniq list
     do j = 1, nbas
        if( abs(dot_product(ax_list(:,i), uniq_ax(1:3,j))) .gt. 0.99_rp ) then
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
  if( count(op .eq. OP_PROP_ROT) .gt. 0 ) then
     max_n_val = maxval( n_int(:), mask=(op .eq. OP_PROP_ROT))
     max_n_loc = maxloc( n_int(:), dim=1, mask=(op .eq. OP_PROP_ROT))
  endif
  !!
  !! special case for D2: n of principal ax is equal to n of other C ax
  !!
  !! D2 has all axes multip 1, any of the C2 axes can be principal,
  !! but D2d has 2 C2 axes which are not prncipal, with multip=1, and principal C2 with multip=3
  if( maxval(multip_ax, mask=(op .eq. OP_PROP_ROT)) .gt. 1 ) then
     !! there are axes with multip > 1
     if( multip_ax( max_n_loc) == 1 ) then
        !! we have chosen the wrong axis, choose again among axes with multip > 1
        max_n_val = maxval( n_int(:), mask=(op .eq. OP_PROP_ROT .and. multip_ax .gt. 1))
        max_n_loc = maxloc( n_int(:), dim=1, mask=(op .eq. OP_PROP_ROT .and. multip_ax .gt. 1))
     endif
  end if

  !! can happen if there are no C operations in PG, for example Cs. Choose the ax of second op
  if( max_n_loc .eq. 0 ) max_n_loc = 2


  ! if(verb) then
  !   write(*,*) 'largest n:', max_n_val, max_n_loc
  !   write(*,*) 'ax:',ax_list(:,max_n_loc)
  ! endif

  !! get how many C ax have this n
  nr_n = count( nint(uniq_ax(4,:)) .eq. max_n_val )
  ! write(*,*) 'nr_n', nr_n

  ! write(*,*) "nbas:",nbas
  ! write(*,*) 'has inversion:',has_inversion
  ! write(*,*) 'has sigma:', has_sigma
  ! write(*,*) 'has cn:', has_cn
  ! write(*,*) 'max_n_val:',max_n_val
  ! write(*,*) 'max_n_loc',max_n_loc
  ! write(*,'(a14,3f9.4)') 'principal ax:',ax_list(:,max_n_loc)

  !! select ax with largest n
  cn_ax = ax_list(:,max_n_loc)

  !! set principal ax for output
  ! prin_ax = cn_ax

  !!====
  !! find all equivalent principal axes:
  !! proper rotation, same n value, same multiplicity
  cn_val = max_n_val
  cn_multip = multip_ax(max_n_loc)
  n_prin_ax = 0
  do i = 1, nbas
     !! find which ax on list have same characteristics as cn
     !! NOTE: this will find duplicates also
     if( n_int(i) .eq. cn_val .and. &
          op(i) .eq. OP_PROP_ROT .and. &
          multip_ax(i) .eq. cn_multip ) then
        ! write(*,"(a,1x, 3f9.4)") "ax is prin:", ax_list(:,i)
        !! store only unique
        isnew = .true.
        do j = 1, n_prin_ax
           if( dot_product(ax_list(:,i), prin_ax(:,j)) .gt. 0.999_rp ) isnew = .false.
        end do
        if( .not. isnew ) cycle
        n_prin_ax = n_prin_ax + 1
        prin_ax(:,n_prin_ax) = ax_list(:,i)
     end if
  end do

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
  dd = 0.0_rp
  n_c2 = 0
  ! write(*,*) '>> checking c2 perp'
  do i = 1, nbas
     if( op(i) .ne. OP_PROP_ROT) cycle
     if( n_int(i) .ne. 2) cycle
     if( i .eq. max_n_loc ) cycle
     ax = ax_list(:,i)
     !! keep dot prod of all
     dot = abs( dot_product( cn_ax, ax))
     !! can be same ax as cn_ax
     if( dot .gt. 0.99_rp) cycle
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
  if( dd .lt. 0.01_rp .and. n_c2 .gt. 0 ) c2_perp = .true.
  ! write(*,*) '>> c2 are perp',c2_perp

  !! check if has sigma_h
  !! horizontal refl plane: axis of plane is parallel to principal, dot( ax_s0, principal ) = 1.0
  ! write(*,*) '>> checking sigma_h'
  do i = 1, nbas
     if( op(i) .ne. OP_IMPROP_ROT) cycle
     if( n_int(i) .ne. 0) cycle
     dot = abs( dot_product(cn_ax, ax_list(:,i)))
     if( dot .gt. 0.99_rp ) has_sigma_h = .true.
     ! write(*,*) i, dot
  end do
  ! write(*,*) '>> sigma_h found:',has_sigma_h



  !! check if has n sigma_v:
  !! vertical plane, contains the principal ax: the norm of cross( ax_s0, principal ) is 1.0
  ! write(*,*) '>> checking sigma_v'
  n_in_plane = 0
  do i = 1, nbas
     if( op(i) .ne. OP_IMPROP_ROT) cycle
     if( n_int(i) .ne. 0) cycle
     !! cross prod( ax, cn_ax ) = 1 for ax in plane
     call cross_prod( ax_list(:,i), cn_ax, cc)
     cross=norm2(cc)
     ! write(*,'(5x,i2,f5.2,3f9.4)') i, cross, ax_list(:,i)
     if( abs( cross ) .gt. 0.99_rp ) n_in_plane = n_in_plane + 1
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
        if( op(i) .ne. OP_PROP_ROT ) cycle
        if( n_int(i) .ne. 2) cycle

        dot = abs( dot_product( ax_list(:,i), cn_ax))
        !! i is perp to cn
        if( dot .lt. 0.01_rp) then
           !! select another
           do j = 1, nbas
              if( j .eq. i ) cycle
              if( op(j) .ne. OP_PROP_ROT ) cycle
              if( n_int(j) .ne. 2) cycle

              dotj = abs( dot_product( ax_list(:,j), cn_ax))
              if( dotj .lt. 0.01_rp) then
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
                    if( dotk .lt. 0.01_rp .and. abs(cross) .gt. 0.99_rp ) then
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
                    if( op(k) .ne. OP_IMPROP_ROT) cycle
                    if( n_int(k) .ne. 0) cycle
                    dotk = abs( dot_product(ax, ax_list(:,k)) )
                    call cross_prod( ax_list(:,k), cn_ax, cc)
                    cross=norm2(cc)
                    ! write(*,*) i,j,k, dotk
                    if( dotk .lt. 0.01_rp .and. cross.gt.0.99_rp ) then
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
                    if( op(k) .ne. OP_IMPROP_ROT) cycle
                    if( n_int(k) .ne. 0) cycle
                    dotk = abs( dot_product(ax, ax_list(:,k)) )
                    call cross_prod( ax_list(:,k), cn_ax, cc)
                    cross=norm2(cc)
                    ! write(*,*) i,j,k, dotk
                    if( dotk .lt. 0.01_rp .and. cross.gt.0.99_rp) then
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
                    if( op(k) .ne. OP_IMPROP_ROT) cycle
                    if( n_int(k) .ne. 0) cycle
                    dotk = abs( dot_product(ax, ax_list(:,k)) )
                    call cross_prod( ax_list(:,k), cn_ax, cc)
                    cross=norm2(cc)
                    ! write(*,*) i,j,k, dotk
                    if( dotk .lt. 0.01_rp .and. cross.gt.0.99_rp) then
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
     if( op(i) .ne. OP_IMPROP_ROT ) cycle
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

  !! Cs and Ci have no prin_ax up to now
  if( pg == "Cs" .or. pg == "Ci" )then
     !! set prin_ax as ax of second Op
     n_prin_ax = 1
     prin_ax(:,1) = ax_list(:,2)
  end if


100 continue

  if( verb) then
     write(*,*) " >> Report from get_pg:"
     write(*,*) 'has inversion:',has_inversion
     write(*,*) 'has sigma:', has_sigma
     write(*,*) 'has cn:', has_cn
     write(*,*) 'max_n_val:',max_n_val
     write(*,'(a,1x,i0)') 'principal ax:', n_prin_ax
     do i = 1, n_prin_ax
        write(*,"(3f9.4)") prin_ax(:,i)
     end do
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
  if( pg_err .eq. 1 ) then
     write(pg,'(a,a1)') trim(pg),'+'
  elseif( pg_err .eq. -1 )then
     write(pg,'(a,a1)') trim(pg),'-'
  elseif( pg_err .eq. -9 ) then
     !! error in check_pg_Nop
     ierr = -1
     write(*,*) "at:",__FILE__," line:",__LINE__
     return
  endif

  ! write(*,*) 'PG is: ',pg

end subroutine sofi_get_pg

!> @brief
!! Analyse the input 3x3 matrix rmat, return the Schoenflies PG
!! notation: "Op n^p", also give axis, angle of the operation.
!!
!! @details
!! The matrices \f$M\f$ corresponding to PG
!! symmetry operations are orthonormal (\f$M^{-1}=M^T\f$), non-symmetric matrices
!! with real elements, and correspond to either pure rotations (symbol \f$C\f$),
!! pure reflections (symbol \f$\sigma\f$), or combinations of both (rotoreflections, symbol \f$S\f$).
!!
!! The matrix properties are:
!! * for C matrices: \f$det=1\f$, \f$\lambda_1 = 1\f$, and \f$\lambda_{2,3}=\exp{(\pm i \alpha)}\f$,
!! * for \f$\sigma\f$ matrices: \f$det=-1\f$, \f$\lambda_1=-1\f$, and \f$\lambda_{2,3}=1\f$.
!! * for S matrices: \f$det=-1\f$, \f$\lambda_1=-1\f$, and \f$\lambda_{2,3}=\exp{(\pm i\alpha)}\f$.
!!
!! The axis always corresponds to the eigenvector of \f$\lambda_1\f$.
!!
!! The approach to characterise the matrices used in this routine is based on value of
!! their determinant, eigenvalues, and eigenvectors.
!!
!! The angle \f$\alpha\f$ in radians is transformed to integers \f$n\f$ and \f$p\f$,
!! such that \f${\alpha}/{2\pi} = {p}/{n}\f$.
!!
!! The use of \f$\sigma\f$ symbol is redundant, since an equivalent notation is a
!! rotoreflection \f$S\f$ with angle \f$\alpha=0\f$. We thus label pure reflections as \f$S~0\f$.
!!
!! The identity E and point-inversion I matrices always have the same form,
!! E=diag(1,1,1) and I=diag(-1,-1,-1). Also, E is equivalent to C with angle zero, and I is equivalent to
!! S with angle 0.5, regardless of the axis.
!!
!! @warning
!!   Axis output from this routine follows the convention:
!!     (1) flip ax such that z>0
!!     (2) if z==0, then flip such that x>0
!!     (3) if x==0, then flip such that y>0.
!!   The value of angle is then relative to this axis convention.
!!
!! @note
!!   The convention for axis oriontation is established due to the reason
!!   that two matrices M and M^T can produce the same result from the
!!   diagonalisation procedure in this routine, meaning the
!!   relative orioentation (positive/negative angle, axis direction)
!!   can be ambiguous.
!!
!! @note
!!   A different method (potentially faster) to get the angles could be based on the trace of matrix,
!!   which should be 3 for E, -3 for I, 1 for s, 2cos(a)+1 for C, 2cos(a)-1 for S.
!!   The axis can be obtained by summing all transformations of a generic vector with all matrices
!!   of the same operation and order (all p from 1 to n, or 2n for S_odd). Multiply by det(M) for S cases.
!!
!! @param[in] rmat(3,3) :: the input matrix
!! @param[out] op(1) :: character for operation E, I, C, S
!! @param[out] n :: the order of operation
!! @param[out] p :: the power for C5^2 kinda things
!! @param[out] ax(3) :: the axis of operation, according to convention
!! @param[out] angle :: angle in units of 1/2pi, e.g. value 0.5 is half the circle
!! @param[out] ierr :: error value, zero on normal execution, negative otherwise
!! @returns op, n, p, ax, angle, ierr
!!
subroutine sofi_analmat( rmat, op, n, p, ax, angle, ierr )
  use ira_precision
  use sofi_tools
  use err_module
  implicit none
  real(rp), dimension(3,3), intent(in) :: rmat
  character(len=1),     intent(out) :: op
  integer(ip),              intent(out) :: n
  integer(ip),              intent(out) :: p
  real(rp), dimension(3),   intent(out) :: ax
  real(rp),                 intent(out) :: angle
  integer(ip),              intent(out) :: ierr

  real(rp), dimension(3,3) :: rdum
  real(rp) :: det, search_eval, diff, diff_old, cosangl
  real(rp), dimension(3) :: eigvals
  integer(ip) :: idx, i, j, nl, pl, gcd
  real(rp) :: mindiff
  ! real(rp) :: flip
  real(rp), dimension(3) :: tmp, tmp2, tmp3

  op = OP_ERROR
  n = 0; p =0


  !! working copy
  rdum(:,:) = rmat(:,:)

  !! determinant
  !! IRA routine
  call determinant3x3( rdum, det )
  if( abs(det) .gt. 1.0_rp + epsilon ) then
     write(*,*) "input matrix:"
     write(*,'(3f9.4)') rdum(1,:)
     write(*,'(3f9.4)') rdum(2,:)
     write(*,'(3f9.4)') rdum(3,:)
     write(*,*) "determinant has value:", det
     write(*,*) "origin at",__FILE__,"line:",__LINE__
     ierr = ERR_DETERMINANT
     return
  end if


  !! diagonalise
  call diag( 3, rdum, eigvals, 1)

  ! write(*,*) 'diagonalised'
  ! do i = 1, 3
  !    write(*,'(3f9.4)') rdum(i,:)
  ! end do

  !! based on value of det, decide what to check
  if( det .gt. 0.5_rp ) then
     !! positive 1.0 determinant, matrix is rotation
     !! there is one eignevalue which is 1.0
     search_eval = 1.0_rp
     ! write(*,*) "matrix is rotation"
     op(1:1)=OP_PROP_ROT
     !!
  elseif( det .lt. -0.5_rp) then
     !! negative 1.0 determinant, matrix is (roto-)inversion.
     !! There is one eignevalue -1.0
     search_eval = -1.0_rp
     ! write(*,*) "matrix is (roto-)inversion"
     op(1:1)=OP_IMPROP_ROT
  endif

  !! find requested eigenvalue (they are not ordered)
  diff_old = 99.9_rp
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
  ! if( cosangl .gt. 1.0 .or. cosangl .lt. -1.0) cosangl = real(nint(cosangl))
  if( abs(cosangl) .gt. (1.0_rp + epsilon) ) then
     !! value cosangl beyond 1.0, or below -1.0, error
     write(*,'(a,1x,3f9.6)') "eigvals:", eigvals
     write(*,'(a,1x,f6.2)') "search_eval:", search_eval
     write(*,'(a,1x,f9.6)') "cosangl has value:",cosangl
     write(*,*) "origin at:",__FILE__,"line:",__LINE__
     ierr = ERR_ACOS_ARG
     return
  elseif( abs(cosangl) .gt. 1.0_rp ) then
     !! value 1.0 within precision, take sign with value 1.0
     cosangl = sign(1.0_rp, cosangl )
  end if

  angle = acos( cosangl ) / (2.0_rp*pi)
  !! the resulting angle here is always positive

  ! write(*,*) cosangl, angle
  ! write(*,*) 'analmat angle:',angle*360

  !! Schoenflies notation: "OP n^p"
  !! primitive way to find n and p. There should be a better way for this...!
  n = 0
  p = 1
  nl = n
  pl = p
  mindiff=99.9_rp
  do j = 1, lim_n_val
     diff = angle*j - nint(angle*j)
     ! write(*,*) j, diff
     if( abs(diff) .lt. mindiff ) then
        nl = j
        pl = nint(angle*j)
        ! exit
        mindiff = abs(diff)
     end if
  end do
  !!
  if( nl .eq. n .and. pl .eq. p ) then
     write(*,'(a,1x,f12.6)') "unable to find n and p for angle:",angle
     write(*,*) "origin at:",__FILE__,"line:",__LINE__
     ierr = ERR_OTHER
     return
  end if

  n = nl
  p = pl
  !! greatest common denominator
  gcd = gcd_rec( p, n)
  ! write(*,*) "gcd:",n, p, gcd
  n = n/gcd
  p = p/gcd


  if( angle .lt. epsilon ) n = 0
  if( p .eq. 0 ) p = 1
  !!
  !! C0 = E :: rotation of angle 0 around any axis is identity
  if( op(1:1)==OP_PROP_ROT .and. angle .lt. epsilon ) op=OP_IDENTITY
  !!
  !! S2 = I :: rotation 0.5 and reflection about any axis is inversion
  if( op(1:1)==OP_IMPROP_ROT .and. abs(angle-0.5_rp) .lt. epsilon ) op=OP_INVERSION

  ! write(*,'(a,x,i3,a3,2f9.4)') 'angle:',n, op, angle, (n-1.0/angle)

  !! detection of possible errors
  if( op(1:1)==OP_PROP_ROT .and. n==0) then
     op=OP_ERROR
     ierr = ERR_OTHER
     write(*,*) "Unknown error in sofi_analmat!"
     write(*,*) "origin at:",__FILE__,"line:",__LINE__
     return
  end if

  !! convention for axis direction:
  !! flip such that z>0
  !! if z==0, then flip such that x>0
  !! if x==0, then flip such that y>0
  call ax_convention( ax )

  !! check if angle is positive/negative according to axis:
  !! generate generic off-axis tmp vector
  tmp = ax + (/ax(3)*0.5_rp, ax(1)*0.2_rp, ax(2)/)
  tmp = tmp/norm2(tmp)
  !! transform tmp by rmat
  tmp2 = matmul(rmat, tmp)
  !! cross product the result with initial
  call cross_prod( tmp, tmp2, tmp3 )
  !! if dot( tmp3, ax ) is negative, angle is negative
  !! if it's zero, the only possibility is angle=0.5, in which case M=M^T
  if( dot_product( tmp3, ax ) < -epsilon ) angle = -angle


  !! put ax of E or I ops to (1, 0, 0 )
  if( op(1:1) == OP_IDENTITY )  ax = (/ 1.0_rp, 0.0_rp, 0.0_rp /)
  if( op(1:1) == OP_INVERSION ) ax = (/ 1.0_rp, 0.0_rp, 0.0_rp /)

  ierr = 0

end subroutine sofi_analmat



!> @details
!! Impose the action of an external magnetic field to the symmetry
!! elements of the PG.
!!
!! The effect of external B field on PG is to filter the list of Ops,
!! only those which satisfy:
!!
!!    det( M ) M B = B
!!
!! are valid, where M is the symmOp matrix.
!!
!! See Eq(2) of:
!! A. Pausch, M Gebele, W. Klopper, J. Chem. Phys. 155, 201101 (2021)
!! https://doi.org/10.1063/5.0069859
!!
!! @param[inout] n_op :: number of symmetry operations
!! @param[inout] op_list(3,3,n_op) :: the list of symmetry operations
!! @param[in] b_field(3) :: direction of the B field
!! @returns n_op, op_list
!!
subroutine sofi_ext_Bfield( n_op, op_list, b_field )
  use ira_precision
  use sofi_tools, only: op_valid_ext_b
  implicit none
  integer(ip),                           intent(inout) :: n_op
  real(rp), dimension(1:3, 1:3, 1:n_op), intent(inout) :: op_list
  real(rp), dimension(3),                intent(in) :: b_field

  integer(ip) :: i, n_op_new
  integer(ip), dimension(n_op) :: which
  real(rp), dimension(3,3) :: rmat
  real(rp), dimension(3,3,n_op) :: op_new

  logical :: valid

  !! array to keep track of which ops survive the filter
  which(:) = 0

  n_op_new = 0
  op_new(:,:,:) = 0.0_rp

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
  op_list(:,:,:) = 0.0_rp
  op_list(:,:,1:n_op) = op_new(:,:,1:n_op)

  if( n_op .lt. 1) then
     !! should not, at least E should satisfy!
     write(*,*) 'heavy error in sofi_ext_Bfield'
     stop
  endif

  return
end subroutine sofi_ext_Bfield


!> @details
!! construct a 3x3 matrix from input data of Op, axis, and angle.
!!
!! A pure rotation matrix \f$C\f$ can be constructed from knowing the angle and
!! axis of rotation, for example by the well-known Euler-Rodrigues rotation formula.
!! See: https://www.wikipedia/rotation_matrix#quaternion
!!
!! A reflection matrix \f$\sigma\f$ can be constructed from knowing the normal vector of
!! the plane of reflection \f$\ket{x}\f$ by \f$\sigma=E - 2\ket{x}\bra{x}\f$, where \f$E\f$ is
!! the identity matrix, and \f$\ket{.}\bra{.}\f$ indicates the outer product.
!!
!! Matrices corresponding to rotoreflections \f$S\f$ can be constructed from combination of pure
!! rotation and pure reflection matrices, \f$S = \sigma C\f$, where \f$\sigma\f$ and \f$C\f$
!! have the same vector as the axis of rotation, and the normal of the plane of reflection.
!!
!! @param[in] op(1) :: character of symmetry operation E, I, C, S
!! @param[in] axis(3) :: the axis
!! @param[in] angle :: the angle in units 1/(2pi), e.g. angle=0.5 means half circle
!! @param[out] matrix(3,3) :: the orthonormal matrix
!! @param[out] ierr :: error value, negative on error, zero otherwise
!!
subroutine sofi_construct_operation( op, axis, angle, matrix, ierr )
  use ira_precision
  use sofi_tools, only: construct_rotation, construct_reflection
  use err_module
  implicit none
  character(len=1),     intent(in) :: op
  real(rp), dimension(3),   intent(in) :: axis
  real(rp),                 intent(in) :: angle
  real(rp), dimension(3,3), intent(out) :: matrix
  integer(ip),              intent(out) :: ierr

  real(rp), dimension(3,3) :: tmp
  real(rp), dimension(3) :: ax
  real(rp) :: an
  real(rp), parameter :: twopi=8.0_rp*atan(1.0_rp)

  ierr = 0
  !! put angle into 2pi units
  an = angle*twopi
  !! normalize axis
  ax = axis/norm2(axis)
  !! initialize
  matrix(:,:) = 0.0_rp
  !!
  select case( trim(adjustl(op)) )
  case( 'E', 'e' )
     !! identity
     matrix(1,1) = 1.0_rp
     matrix(2,2) = 1.0_rp
     matrix(3,3) = 1.0_rp
  case( 'I', 'i' )
     !! inversion
     matrix(1,1) = -1.0_rp
     matrix(2,2) = -1.0_rp
     matrix(3,3) = -1.0_rp
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
     write(*,*) "origin at:",__FILE__,"line:",__LINE__
     ierr = ERR_UNKNOWN_OP
     return
  end select

end subroutine sofi_construct_operation

subroutine sofi_unique_ax_angle( n_mat, mat_list, op_out, ax_out, angle_out, ierr )
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

  use ira_precision
  use sofi_tools
  use err_module
  implicit none
  integer(ip),                    intent(in) :: n_mat
  real(rp), dimension(3,3,n_mat), intent(in) :: mat_list
  character(len=1), dimension(n_mat), intent(out) :: op_out
  real(rp), dimension(3,n_mat),           intent(out) :: ax_out
  real(rp), dimension(n_mat),             intent(out) :: angle_out
  integer(ip),                            intent(out) :: ierr

  integer(ip) :: n, i, p, j
  real(rp) :: angle, dotp, angle_diff, angle_sum, dist, dist_neg
  character(len=1), dimension(n_mat) :: op_list
  real(rp), dimension(3) :: ax
  real(rp), dimension(3,3) :: rmat
  real(rp) :: dotp_equal

  ierr = 0

  !! threshold value for dot product between two vectors which should be equal
  dotp_equal = 0.9999_rp

  !! get all axes
  do i = 1, n_mat
     rmat = mat_list(:,:,i)
     call sofi_analmat( rmat, op_list(i), n, p, ax, angle, ierr )
     if( ierr /= 0 ) return
     op_out(i) = op_list(i)
     ax_out(:,i) = ax
     angle_out(i) = angle
     ! write(*,'(i4,x,a2,x,i0,"^",i0,x,3f9.4,4x,f9.4)') i, op_out(i),n,p, ax_out(:,i), angle_out(i)
  end do

  !! find operations that are ambiguous:
  !!   - same ax, same angle
  !!   - opposite ax, opposite angle
  do i = 1, n_mat
     do j = i+1, n_mat
        if( op_list(i) .ne. op_list(j) ) cycle
        !! axes
        dotp=dot_product( ax_out(:,i), ax_out(:,j))
        !! angles
        angle_diff = angle_out(i) - angle_out(j)
        angle_sum = angle_out(i) + angle_out(j)
        !!
        !! same ax, same angle
        if( dotp .gt. dotp_equal .and. &
             abs(angle_diff) .lt. 1e-3_rp ) then
           !! ops are ambiguous
           write(*,*) i, j
           write(*,*)"dotp",dotp, "angle_diff",abs(angle_diff)
           write(*,'(3f9.5,2x,f7.4)') ax_out(:,i), angle_out(i)
           write(*,'(3f9.5,2x,f7.4)') ax_out(:,j), angle_out(j)
           ierr = -1
        end if
        !!
        !! opposite ax, equal or opposite angle
        if( dotp .lt. -dotp_equal  )then
           if( abs( angle_diff ) .lt. 1e-3_rp .or. &
                abs(angle_sum) .lt. 1e-2_rp ) then
              !! ops are ambiguous
              write(*,*) i, j, angle_diff
              write(*,*)"dotp",dotp, "angle_diff",abs(angle_diff), "angle_sum",abs(angle_sum)
              write(*,'(3f9.5,2x,f7.4)') ax_out(:,i), angle_out(i)
              write(*,'(3f9.5,2x,f7.4)') ax_out(:,j), angle_out(j)
              ierr = -1
           endif
        end if
     end do
  end do


  !! put ax of E or I ops to (1, 0, 0 )
  do i = 1, n_mat
     if( op_out(i) == OP_IDENTITY )  ax_out(:,i) = (/ 1.0_rp, 0.0_rp, 0.0_rp /)
     if( op_out(i) == OP_INVERSION ) ax_out(:,i) = (/ 1.0_rp, 0.0_rp, 0.0_rp /)
  end do



  !! output final list
  ! write(*,*) '----- final'
  ! do i = 1, n_mat
  !    write(*,'(i4,x,a,x,3f9.4,4x,f9.4)') i, op_out(i), ax_out(:,i), angle_out(i)
  ! end do


  ! ierr=ERR_OTHER

  !! reconstruct each from op, angle, ax
  !! reconstruct + angle

  !! reconstruct - angle
end subroutine sofi_unique_ax_angle


subroutine sofi_mat_combos( n_in, mat_in, n_out, mat_out )
  !! generate combos of matrices without a structure
  use ira_precision
  use sofi_tools, only: nmax, matrix_distance, m_thr
  implicit none
  integer(ip), intent(in) :: n_in
  real(rp), dimension(3,3,n_in), intent(in) :: mat_in
  integer(ip), intent(out) :: n_out
  real(rp), dimension(3,3,nmax), intent(out) :: mat_out

  integer(ip) :: i, m, ii, j, nn, k
  real(rp), dimension(3,3) :: rmat
  logical :: is_new
  real(rp) :: dd

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

              if( dd .lt. m_thr ) then
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

!> @cond SKIP
!! @param[in]
subroutine state_update( n_old, n_op, op_list, has_sigma, has_inversion, n_u_c3, u_c3, terminate )
  !! update state flags, terminate if condition for Ih is found.
  !! condition for Ih: has sigma or inversion, and number of unique c3 axes .ge. 3
  !!
  !! NOTE: The variables which store the state need to be allocated and initialised in the caller!
  use ira_precision
  use sofi_tools, only: OP_PROP_ROT, OP_IMPROP_ROT, OP_INVERSION
  implicit none
  integer(ip), intent(in)    :: n_old   !! size of previous call (increments are not always +1)
  integer(ip), intent(in)    :: n_op
  real(rp),    intent(in)    :: op_list(3,3,n_op)
  logical, intent(inout) :: has_sigma          !! state var: sigma is present, init to false
  logical, intent(inout) :: has_inversion      !! state var: inversion is present, init to false
  integer(ip), intent(inout) :: n_u_c3             !! state var: number of unique c3 axes, init to zero
  real(rp),    intent(inout) :: u_c3(3,5)          !! state var: unique c3 axes, init to zero
  logical, intent(out)   :: terminate

  !!
  integer(ip) :: i, iax, ierr
  integer(ip) :: n_new
  real(rp), dimension(3,3) :: mat
  real(rp), dimension(3) :: ax
  character(len=1) :: op
  integer(ip) :: n, p
  real(rp) :: angle, dotp
  logical :: unique

  terminate = .false.

  !! check only the new operations
  n_new = n_op - n_old

  do i = 1, n_new
     mat = op_list(:,:,n_old + i)
     call sofi_analmat( mat, op, n, p, ax, angle, ierr )
     if( ierr /= 0 ) then
        write(*,*) "got err from:",__FILE__,__LINE__
        write(*,*) "ierr value:", ierr
        return
     end if
     !!
     !! flip flags
     if( .not. has_sigma) then
        if( op == OP_IMPROP_ROT .and. n==0 ) has_sigma = .true.
     end if
     if( .not. has_inversion ) then
        if( op == OP_INVERSION ) has_inversion = .true.
     end if

     !! if C3, check the unique axes
     c3ax: if( op == OP_PROP_ROT .and. n .ge. 3) then
        !!
        if( n_u_c3 .ge. 5 ) exit c3ax
        unique = .true.
        do iax = 1, n_u_c3
           dotp = dot_product( ax, u_c3(:,iax))
           if( abs(dotp) .gt. 0.99_rp ) unique = .false.
        end do
        !!
        if( unique ) then
           !! add ax to list of unique c3 axes
           n_u_c3 = n_u_c3 + 1
           u_c3(:, n_u_c3 ) = ax
        end if
     end if c3ax

     !! check condition for Ih:
     if( n_u_c3 .ge. 5 ) then
        if( has_sigma .or. has_inversion ) terminate = .true.
     end if

  end do


end subroutine state_update
!> @endcond


!> @details
!! extract the error message from ira/sofi corresponding to error value ierr
!!
!! @param[in] ierr :: integer error value
!! @param[out] msg(128) :: error message
!!
subroutine sofi_get_err_msg( ierr, msg )
  use ira_precision
  use err_module
  implicit none
  integer(ip), intent(in) :: ierr
  character(len=128), intent(out) :: msg
  character(:), allocatable :: str

  !! cant output allocatable ... write to fixed length string
  str = get_err_msg( ierr )
  if( len(str) .le. len(msg) ) then
     msg = str
  else
     write(*,*) "WARNING: error str longer than given msg variable!"
     write(*,*) len(msg), len(str)
     write(*,*) "in: ",__FILE__, __LINE__
     msg = str(1:len(msg))
  end if
end subroutine sofi_get_err_msg


!> @details
!! check if the atomic positions in coords form a linear structure.
!!
!! @param[in] integer(ip) :: nat --> number of atoms
!! @param[in] real(3,nat) :: coords --> atomic positions
!! @param[out] logical :: collinear --> .true. when structure is collinear
!! @param[out] real(3) :: ax --> axis of collinearity
!!
subroutine sofi_check_collinear( nat, coords, collinear, ax_o )
  !! Method: vectors connecting pairs of coords must be collinear
  use ira_precision
  use sofi_tools, only: collinearity_thr, ax_convention
  implicit none
  integer(ip),                intent(in) :: nat
  real(rp), dimension(3,nat), intent(in) :: coords
  logical,                intent(out) :: collinear
  real(rp), dimension(3),     intent(out) :: ax_o

  !! local vars
  real(rp) :: dotp, ax(3), vec(3)
  integer(ip) :: i

  ax_o(:) = 0.0_rp

  collinear = .true.
  if( nat == 1 ) return

  !! choose ax as normalized vector from 1 to 2
  ax_o = coords(:,2) - coords(:,1)
  !! flip if needed
  call ax_convention(ax_o)
  ax_o = ax_o/norm2(ax_o)
  !! keep this as reference ax
  ax = ax_o

  !! nat==2 is collinear for sure
  if( nat == 2 ) return

  ! ax_o = ax_o + ax
  ! write(*,*) 0, ax
  check: do i = 3, nat
     !! dot_product of nomralized vec coords(:,i)-coords(:,1) with ax should be 1
     vec = coords(:,i) - coords(:,i-1)
     vec = vec/norm2(vec)
     !!
     dotp = dot_product( ax, vec )
     ax_o = ax_o + sign(1.0_rp, dotp)*vec
     ! write(*,*) i, vec, dotp
     if( abs(dotp) < collinearity_thr ) then
        collinear = .false.
        exit check
     end if
  end do check

  ax_o = ax_o / norm2(ax_o)

  ! write(*,*) "collinear:",collinear
  ! write(*,*) "ax_o", ax_o
end subroutine sofi_check_collinear
