!!
!! Copyright (C) 2021, MAMMASMIAS Consortium
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


!> @cond SKIP
  module ira_pbc

    !! routines for computing a vector in periodic boundary condtions
    use ira_precision
    implicit none
    private
    public :: pbc_vec

    contains

      subroutine pbc_vec( vec, lat )
        !! apply pbc of lattice 'lat' to a vector 'vec'
        implicit none
        real(rp), dimension(3), intent(inout) :: vec
        real(rp), dimension(3,3), intent(in) :: lat

        call cart_to_crist( vec, lat )
        call periodic( vec )
        call crist_to_cart( vec, lat )

      end subroutine pbc_vec

      subroutine periodic(c)
        !--------------------------------
        ! periodic boundary condition, for 3 dimensional vector input in crist coords.
        !--------------------------------
        implicit none
        real(rp), dimension(3),intent(inout) :: c
        integer(ip) :: i

        do i = 1, 3
          if( c(i) .lt. -0.5_rp ) c(i) = c(i) + 1.0_rp
          if( c(i) .ge. 0.5_rp ) c(i) = c(i) - 1.0_rp
        end do

      end subroutine periodic


      subroutine cart_to_crist(xpp,ct)
        !!----------------------------
        !! cartesian to crystallographic coordinates transform, in 3-dimension
        !! v_crist = B^-1 * R_cart; where B is the matrix formed by unit cell vectors
        !! This routine does the transpose of B implicitly
        !! --------
        !! xpp(3)      ==> input vector of position in cartesian
        !! ct(3,3)     ==> conversion matrix, vectors of the Bravais lattice in rows
        !!
        !!      ct = a1 a2 a3
        !!           b1 b2 b3
        !!           c1 c2 c3
        !!----------------------------
        !! bt(3,3) ==> inverse matrix of ct, used locally
        !! xc(3)   ==> copy of xpp, used locally
        !! detct   ==> determinant of ct, used locally
        !!
        implicit none
        real(rp), dimension(3),   intent(inout) :: xpp
        real(rp), dimension(3,3), intent(in)    :: ct

        real(rp),dimension(3) :: xc
        real(rp) :: detct
        real(rp), dimension(3,3) :: bt

        ! -----------------------------------------------
        !  inverse matrix of ct(:,:)
        !------------------------------------------------
        detct=ct(1,1)*ct(2,2)*ct(3,3)+&
              ct(1,2)*ct(2,3)*ct(3,1)+&
              ct(2,1)*ct(3,2)*ct(1,3)&
              -ct(1,3)*ct(2,2)*ct(3,1)&
              -ct(3,2)*ct(2,3)*ct(1,1)&
              -ct(1,2)*ct(2,1)*ct(3,3)

        bt(1,1)= ct(2,2)*ct(3,3)-ct(2,3)*ct(3,2)
        bt(1,2)=-(ct(1,2)*ct(3,3)-ct(1,3)*ct(3,2))
        bt(1,3)= ct(1,2)*ct(2,3)-ct(1,3)*ct(2,2)
        bt(2,1)=-(ct(2,1)*ct(3,3)-ct(2,3)*ct(3,1))
        bt(2,2)= ct(1,1)*ct(3,3)-ct(3,1)*ct(1,3)
        bt(2,3)=-(ct(1,1)*ct(2,3)-ct(1,3)*ct(2,1))
        bt(3,1)= ct(2,1)*ct(3,2)-ct(2,2)*ct(3,1)
        bt(3,2)=-(ct(1,1)*ct(3,2)-ct(1,2)*ct(3,1))
        bt(3,3)= ct(1,1)*ct(2,2)-ct(2,1)*ct(1,2)
        !------------------------------------------------

        xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))/detct
        xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))/detct
        xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))/detct

        xpp(:) = xc(:)

      end subroutine cart_to_crist


      subroutine crist_to_cart(xpp,bt)
        !!--------------------------------
        !! crystallographic to cartesian transformation in 3-dimensions
        !! R_cart = B * v_crist; where B is the matrix formed by cell vectors vertically
        !! This routine does the transpose implicitly!
        !! -----------
        !! xpp(3)    ==> input vector in crystallographic, output vector in cartesian
        !! bt(3,3)   ==> input conversion matrix, vectors of the Bravais lattice in rows
        !!
        !!         bt = a1 a2 a3
        !!              b1 b2 b3
        !!              c1 c2 c3
        !!-----
        !! xc(3)   ==> local vector
        !!
        implicit none

        real(rp), dimension(3),   intent(inout) :: xpp

        real(rp), dimension(3,3), intent(in)    :: bt
        real(rp), dimension(3) :: xc


        xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))
        xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))
        xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))

        xpp(:) = xc(:)

      end subroutine crist_to_cart

  end module
!> @endcond



  !> @details
  !! Linear Assignment Problem (LAP) algorithm:
  !! Constrained Shortest Distance Assignment (CShDA).
  !!
  !! Assign atoms of configuration 1 to atoms of configuration 2 based on
  !! the CShDA algorithm.
  !!
  !! @note
  !!  WARNING nat1 must be lower or equal nat2!
  !!
  !!================================================
  !! @param[in]    nat1    -> number of atoms in conf 1;
  !! @param[in]    typ1    -> atomic types in conf 1;
  !! @param[in]    coords1 -> coordinates of conf 1;
  !! @param[in]    nat2    -> number of atoms in conf 2;
  !! @param[in]    typ2    -> atomic types in conf 2;
  !! @param[in]    coords2 -> coordinates of conf 2;
  !! @param[in]    some_threshold -> threshold for the Hausdorff distance, used for early exit;
  !! @param[out]   found   -> list of assigned atoms of conf 2 to conf 1:
  !!                          e.g. found(3) = 9 means atom 3 from conf 1 is assigned
  !!                          to atom 9 in conf 2;
  !! @param[out]   dists   -> distances from atom i in conf 1 to atom found(i) in conf 2;
  !!
  subroutine cshda( nat1, typ1, coords1, &
                    nat2, typ2, coords2, &
                    some_threshold, found, dists )
    use ira_precision
    implicit none
    integer(ip),                  intent(in) :: nat1
    integer(ip), dimension(nat1), intent(in) :: typ1
    real(rp), dimension(3,nat1),  intent(in) :: coords1
    integer(ip),                  intent(in) :: nat2
    integer(ip), dimension(nat2), intent(in) :: typ2
    real(rp), dimension(3,nat2),  intent(in) :: coords2
    real(rp),                     intent(in) :: some_threshold
    integer(ip), dimension(nat2), intent(out) :: found
    real(rp), dimension(nat2),    intent(out) :: dists
    !!
    !! local
    !!
    real(rp), dimension(3) :: rij
    real(rp), dimension(3) :: ci, cj
    real(rp) :: dx, dy, dz, th, m_th, th2
    real(rp), dimension(nat2,nat1) :: chkmat
    logical, dimension(nat1) :: lsearch
    integer(ip), dimension(nat2) :: assigned
    integer(ip) :: i, j, k, ti, tj
    integer(ip) :: idx_old
    integer(ip) :: n_count, nmax
    real(rp) :: dist, dist_old, dmin
    integer(ip), dimension(nat1) :: tmpmin


    !! init output
    dists = 999.9_rp
    found = 0

    !!
    !! set up distance matrix, compute elements as dist^2, do sqrt at the end
    !!
    !! ATTENTION: the double loop needs to go from 1 to nat, for i and j
    !!    because there are permutations in the set, so the distance matrix is
    !!    not symmetric!!!
    !!
    th = some_threshold
    th2 = th*th
    m_th = -th
    do i = 1, nat1
       ti = typ1(i)
       ci = coords1(:,i)
       dmin = huge( dmin )
       do j = 1, nat2
          !!
          !! if the atoms are not of same typ, set large value for distance:
          !! like this they will not be found assigned
          !!
          if( ti .eq. typ2(j) ) then
             !!
             cj = coords2(:,j)
             !!
             dx = ci(1) - cj(1)
             dy = ci(2) - cj(2)
             dz = ci(3) - cj(3)
             !!
             !! dist squared
             dist = dx*dx + dy*dy + dz*dz
             !!
             chkmat(j,i) = dist
             !!
             !! keep for minimum row
             ! dmin = min( dmin, dist )
             if( dist .lt. dmin ) then
                dmin = dist
                tmpmin(i) = j
             end if
          else
             chkmat(j,i) = 995.0_rp
          end if
          !!
       end do
       !!
       !! Early exit:
       !! if any row chkmat(i,:) has all values above some_threshold,
       !! then there is no way that Hausdorff distance be lower than some_threshold.
       !! This criterion is used for early return of cshda.
       if( dmin .gt. th2 ) then
          return
       endif
       !!
    end do

    ! write(*,"(3x)",advance="no")
    ! do i = 1, nat2
    !    write(*,"(i4,1x)", advance="no") i
    ! end do
    ! write(*,*)
    ! do i = 1, nat1
    !    write(*,"(i2,1x,*(f4.2,:,1x))")  i, chkmat(:,i)
    ! end do


    !! set up the queue of searches
    lsearch(:) = .true.

    assigned(:) = 0

    !! set first search index
    i = 1

    n_count = 1
    nmax = nat1*nat2
    !! in worst case do n searches on each of the n sites -> nmax = n**2
    do n_count = 1, nmax
       !!
       !! set index of next search
       !!
       i = findloc( lsearch(:), .true., 1)
       if( i .eq. 0 ) exit
       ! write(*,*) "next i",i
       !!
       !! set next search on this index to .false.
       !!
       lsearch(i) = .false.
       !!
       !!
       !! find minimum distance and its index
       !!
       j = tmpmin(i)
       dist = chkmat(j,i)
       ! write(*,*) "tmpmin j", j, dist
       !!
       !!
       !! check the if we already have this j
       !!
       ! if( assigned(j) .gt. 0 ) then
       if( any(found .eq. j) ) then
          !!
          !!
          !! find the old index where its used, and the old distance
          !!
          ! idx_old = assigned(j)
          ! idx_old = minloc( abs(found - j), 1)
          idx_old = findloc( found, j, 1)
          dist_old = dists(idx_old)
          ! write(*,*) "j already assigned at", idx_old, dist_old
          ! write(*,"(*(f4.2,:,1x))") dists
          !!
          if( dist_old .lt. dist ) then
             !!
             !!
             !! if the previous found is closer, set the current distance
             !! to smth big, so its not found ever again!
             !!
             chkmat(j,i) = 999.0_rp
             tmpmin(i) = minloc(chkmat(:,i),1)
             !!
             !!
             !! and the current index should be searched again
             !!
             lsearch(i) = .true.
             !!
          else
             !!
             !!
             !! if the previous found is larger then the new, the old idx should
             !! be searched again and the same distance should not be found!
             !!
             chkmat(j, idx_old) = 999.0_rp
             lsearch( idx_old ) = .true.
             tmpmin(idx_old) = minloc( chkmat(:,idx_old), 1)
             assigned(idx_old) = 0
             !!
          endif
       endif
       !!
       !! set found data
       !!
       found(i) = j
       dists(i) = dist
       assigned(j) = i
       ! write(*,*) "found now"
       ! write(*,"(10i3)") found

       !!
    end do


    !! do sqrt of dists
    do i = 1, nat1
       dists(i) = sqrt(dists(i))
    end do


    !! for equal sizes of structures we should be done
    if( nat1 .eq. nat2 ) return

    !! if any found maps to zero, cshda has exited
    if( any(found(1:nat1) .eq. 0) ) return

    ! write(*,*) "found now"
    ! write(*,"(10i3)") found

    !! find indices of conf2 not represented in found,
    !! and put them at end
    k = nat1
    do i = 1, nat2
       !! if this i is already found, do nothing
       if( any(i .eq. found(:) ) ) cycle
       !! add this i to last spot
       ! write(*,*) "adding ",i,"to idx k=", k
       k = k + 1
       found(k) = i
       ! write(*,*) "found now"
       ! write(*,"(10i3)") found

    end do



  end subroutine cshda


  !> @details
  !! Linear Assignment Problem (LAP) algorithm:
  !! Constrained Shortest Distance Assignment (CShDA).
  !!
  !! Assign atoms of configuration 1 to atoms of configuration 2 based on
  !! the CShDA algorithm. The configuration 2 is periodic with given lattice.
  !!
  !!
  !! @note
  !!  WARNING nat1 must be lower or equal nat2!
  !!
  !!================================================
  !! @param[in]  nat1    -> number of atoms in conf 1;
  !! @param[in]  typ1    -> atomic types in conf 1;
  !! @param[in]  coords1 -> coordinates of conf 1;
  !! @param[in]  nat2    -> number of atoms in conf 2;
  !! @param[in]  typ2    -> atomic types in conf 2;
  !! @param[in]  coords2 -> coordinates of conf 2;
  !! @param[in]  lat2    -> lattice vectors of conf 2 in rows;
  !! @param[in]  some_thr -> threshold for hd;
  !! @param[out]  found   -> list of paired atoms of conf 2 to conf 1:
  !!                         e.g. found(3) = 9 means atom 3 from conf 1 is paired
  !!                         to atom 9 in conf 2;
  !! @param[out]  dists   -> distances from atom i in conf 1 to atom found(i) in conf 2;
  !!
  subroutine cshda_pbc( nat1, typ1, coords1, &
                   nat2, typ2, coords2, lat2, &
                   some_thr, found, dists )

    use ira_precision
    use ira_pbc, only: pbc_vec
    implicit none
    integer(ip),                  intent(in) :: nat1
    integer(ip), dimension(nat1), intent(in) :: typ1
    real(rp), dimension(3,nat1),  intent(in) :: coords1
    integer(ip),                  intent(in) :: nat2
    integer(ip), dimension(nat2), intent(in) :: typ2
    real(rp), dimension(3,nat2),  intent(in) :: coords2
    real(rp), dimension(3,3), intent(in) :: lat2
    real(rp),                     intent(in) :: some_thr
    integer(ip), dimension(nat2), intent(out) :: found
    real(rp), dimension(nat2),    intent(out) :: dists
    !!
    !! local
    !!
    real(rp), dimension(3) :: rij, rj
    real(rp), dimension(nat1,nat2) :: chkmat
    integer(ip), dimension(nat1) :: search
    integer(ip) :: i, j, k
    integer(ip) :: idx_old
    integer(ip) :: n_count
    real(rp) :: dist, dist_old

    dists(:) = 999.9_rp
    !! set found to zero
    found(:) = 0
    !!
    !! set up distance matrix
    chkmat(:,:) = 0.0_rp
    !!
    !! ATTENTION: the double loop needs to go from 1 to nat, for i and j
    !!    because there are permutations in the set, so the distance matrix is
    !!    not symmetric!!!
    !!
    !! compute distances from struc1 to struc2, where struc2 is in pbc
    !!
    do i = 1, nat1
       do j = 1, nat2
          !!
          ! rj(:) = coords2(:,j)
          rij = coords1(:,i) - coords2(:,j)
          call pbc_vec( rij, lat2 )
          dist = sqrt( dot_product(rij, rij) )
          !!
          !! if the atoms are not of same typ, set some large distance:
          !! like this they will not be found paired
          !!
          if( typ1(i) .ne. typ2(j) ) dist = 99990.0_rp
          !!
          chkmat(i,j) = dist
          !!
       end do
       !!
       !! Early return method:
       !! if all values in the row chkmat(i,:) are above some_threshold, there
       !! is no way that final dH could be below that threshold
       if( minval(chkmat(i,:)) .gt. some_thr ) then
          return
       endif
       !!
    end do

    !! set up the queue of searches
    search(:) = 1

    !! set first search
    i = 1
    j = minloc( chkmat(i,:), 1 )

    n_count = 1
    do while( search(i) .gt. 0 )
       !!
       !! return on huge number of searches
       !! ( in worst case do n searches on each of the n sites )
       !!
       if( n_count .gt. nat1**2) then
          found(i) = 0
          dists(i) = 999.9_rp
          write(*,*) " PROBLEM in cshda_pbc: huge number of searches"
          return
       endif
       !!
       !!
       !! set next search on this index to 0
       !!
       search(i) = 0
       !!
       !!
       !! find minimum distance and its index
       !!
       j = minloc( chkmat(i,:), 1)
       dist = chkmat(i,j)
       !!
       !!
       !! check the found indices if we already have this j
       !!
       if( any(found .eq. j) ) then
          !!
          !!
          !! find the old index where its used, and the old distance
          !!
          idx_old = minloc( abs( found - j) , 1)
          dist_old = minval( chkmat(idx_old,:), 1)
          !!
          if( dist_old .lt. dist ) then
             !!
             !!
             !! if the previous found is closer, set the current distance
             !! to smth big, so its not found ever again!
             !!
             chkmat(i,j) = 999.9_rp
             !!
             !!
             !! and the current index should be searched again
             !!
             search(i) = 1
             !!
          else
             !!
             !!
             !! if the previous found is larger then the new, the old idx should
             !! be searched again and the same distance should not be found!
             !!
             chkmat(idx_old, j) = 999.9_rp
             search( idx_old ) = 1
             !!
          endif
       endif
       !!
       !!
       !! set found data
       !!
       found(i) = j
       dists(i) = dist
       !!
       !!
       !! set index of next search
       !!
       i = maxloc( search(:), 1 )
       !!
       n_count = n_count + 1
       !!
    end do

    !! find indices of conf2 not represented in found
    k = nat1
    do i = 1, nat2
       !! if this i is already found, do nothing
       if( any(i .eq. found(:) ) ) cycle
       !! add this i to last spot
       k = k + 1
       found(k) = i
    end do

  end subroutine cshda_pbc


