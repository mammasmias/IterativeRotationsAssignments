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


  subroutine cshda( nat1, typ1, coords1, &
                   nat2, typ2, coords2, &
                   found, dists )

    !> @detail
    !! Linear Assignment Problem (LAP) algorithm:
    !! Constrained Shortest Distance Assignment (CShDA).
    !!
    !! Assign atoms of configuration 1 to atoms of configuration 2 based on
    !! the CShDA algorithm.
    !!
    !! WARNING nat1 must be lower or equal nat2!
    !!
    !!================================================
    !! nat1    -> number of atoms in conf 1;
    !! typ1    -> atomic types in conf 1;
    !! coords1 -> coordinates of conf 1;
    !! nat2    -> number of atoms in conf 2;
    !! typ2    -> atomic types in conf 2;
    !! coords2 -> coordinates of conf 2;
    !! found   -> list of assigned atoms of conf 2 to conf 1:
    !!            e.g. found(3) = 9 means atom 3 from conf 1 is assigned
    !!            to atom 9 in conf 2;
    !! dists   -> distances from atom i in conf 1 to atom found(i) in conf 2;
    !!
    implicit none
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1
    real, dimension(3,nat1),  intent(in) :: coords1
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2
    real, dimension(3,nat2),  intent(in) :: coords2
    integer, dimension(nat2), intent(out) :: found
    real, dimension(nat2),    intent(out) :: dists
    !!
    !! local
    !!
    real, dimension(3) :: rij
    real, dimension(nat1,nat2) :: chkmat
    integer, dimension(nat1) :: search
    integer :: i, j, k
    integer :: idx_old
    integer :: n_count
    real :: dist, dist_old

    found(:) = 0
    dists(:) = 999.9
    !!
    !! set up distance matrix
    chkmat(:,:) = 0.0
    !!
    !! ATTENTION: the double loop needs to go from 1 to nat, for i and j
    !!    because there are permutations in the set, so the distance matrix is
    !!    not symmetric!!!
    !!
    do i = 1, nat1
       do j = 1, nat2
          !!
          rij = coords1(:,i) - coords2(:,j)
          dist = sqrt( dot_product(rij, rij) )
          !!
          !! if the atoms are not of same typ, set some large distance:
          !! like this they will not be found assigned
          !!
          if( typ1(i) .ne. typ2(j) ) dist = 990.0
          !!
          chkmat(i,j) = dist
          !!
       end do
       !!
       !! Early exit idea:
       !! if any row chkmat(i,:) has all values above some_threshold,
       !! then there is no way that Hausdorff distance be lower than some_threshold.
       !! Can use this criterion for early return.
       !!
    end do

    !! set up the queue of searches
    search(:) = 1

    !! set found to zero
    found(:) = 0

    !! set first search
    i = 1
    j = minloc( chkmat(i,:), 1 )

    n_count = 1
    do while( search(i) .gt. 0 )
       !!
       !! return on huge number of searches
       !! ( in worst case do n searches on each of the n sites -> n**2 )
       !!
       if( n_count .gt. nat1*nat2) then
          found(i) = 0
          dists(i) = 999.9
          write(*,*) " PROBLEM in dist_set_reg: huge number of searches"
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
             chkmat(i,j) = 999.9
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
             chkmat(idx_old, j) = 999.9
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
       !! early exit idea:
       !!  if any dist is above some_threshold, then Hausdorff cannot be below it.
       !!
       !!
       !! set index of next search
       !!
       i = maxloc( search(:), 1 )
       !!
       n_count = n_count + 1
       !!
    end do

    !! find indices of conf2 not represented in found,
    !! and put them at end
    k = nat1
    do i = 1, nat2
       !! if this i is already found, do nothing
       if( any(i .eq. found(:) ) ) cycle
       !! add this i to last spot
       k = k + 1
       found(k) = i
    end do

  end subroutine cshda


  subroutine cshda_pbc( nat1, typ1, coords1, &
                   nat2, typ2, coords2, lat2, &
                   found, dists )

    !> @detail
    !! Linear Assignment Problem (LAP) algorithm:
    !! Constrained Shortest Distance Assignment (CShDA).
    !!
    !! Assign atoms of configuration 1 to atoms of configuration 2 based on
    !! the CShDA algorithm. The configuration 2 is periodic with given lattice.
    !!
    !!
    !! WARNING nat1 must be lower or equal nat2!
    !!
    !!================================================
    !! nat1    -> number of atoms in conf 1;
    !! typ1    -> atomic types in conf 1;
    !! coords1 -> coordinates of conf 1;
    !! nat2    -> number of atoms in conf 2;
    !! typ2    -> atomic types in conf 2;
    !! coords2 -> coordinates of conf 2;
    !! lat2    -> lattice vectors of conf 2;
    !! found   -> list of paired atoms of conf 2 to conf 1:
    !!            e.g. found(3) = 9 means atom 3 from conf 1 is paired
    !!            to atom 9 in conf 2;
    !! dists   -> distances from atom i in conf 1 to atom found(i) in conf 2;
    !!
    implicit none
    integer,                  intent(in) :: nat1
    integer, dimension(nat1), intent(in) :: typ1
    real, dimension(3,nat1),  intent(in) :: coords1
    integer,                  intent(in) :: nat2
    integer, dimension(nat2), intent(in) :: typ2
    real, dimension(3,nat2),  intent(in) :: coords2
    real, dimension(3,3), intent(in) :: lat2
    integer, dimension(nat2), intent(out) :: found
    real, dimension(nat2),    intent(out) :: dists
    !!
    !! local
    !!
    real, dimension(3) :: rij, rj
    real, dimension(nat1,nat2) :: chkmat
    integer, dimension(nat1) :: search
    integer :: i, j, k
    integer :: idx_old
    integer :: n_count
    real :: dist, dist_old

    !!
    !! set up distance matrix
    chkmat(:,:) = 0.0
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
          rj(:) = coords2(:,j)
          call cart_to_crist(rj,lat2)
          call periodic(rj)
          call crist_to_cart(rj,lat2)
          rij = coords1(:,i) - rj
          dist = sqrt( dot_product(rij, rij) )
          !!
          !! if the atoms are not of same typ, set some large distance:
          !! like this they will not be found paired
          !!
          if( typ1(i) .ne. typ2(j) ) dist = 99990.0
          !!
          chkmat(i,j) = dist
          !!
       end do
    end do

    !! set up the queue of searches
    search(:) = 1

    !! set found to zero
    found(:) = 0

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
          dists(i) = 999.9
          write(*,*) " PROBLEM in dist_set_reg: huge number of searches"
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
             chkmat(i,j) = 999.9
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
             chkmat(idx_old, j) = 999.9
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


  subroutine periodic(c)
    !--------------------------------
    ! periodic boundary condition, for 3 dimensional vector input in crist coords.
    !--------------------------------
    implicit none
    real, dimension(3),intent(inout) :: c
    integer :: i

    do i = 1, 3
       if( c(i) .lt. -0.5 ) c(i) = c(i) + 1.0
       if( c(i) .ge. 0.5 ) c(i) = c(i) - 1.0
    end do

  end subroutine periodic


  subroutine cart_to_crist(xpp,ct)
    !> @detail
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
    real, dimension(3),   intent(inout) :: xpp
    real, dimension(3,3), intent(in)    :: ct

    real,dimension(3) :: xc
    real :: detct
    real, dimension(3,3) :: bt

       bt(:,:)=0.0
       xc(:) = 0.0

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

       return
  end subroutine cart_to_crist


  subroutine crist_to_cart(xpp,bt)
    !> @detail
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
   real, dimension(3),   intent(inout) :: xpp
   real, dimension(3,3), intent(in)    :: bt
   real, dimension(3) :: xc

       xc(:) = 0.0

       xc(1) = (xpp(1)*bt(1,1)+xpp(2)*bt(2,1)+xpp(3)*bt(3,1))
       xc(2) = (xpp(1)*bt(1,2)+xpp(2)*bt(2,2)+xpp(3)*bt(3,2))
       xc(3) = (xpp(1)*bt(1,3)+xpp(2)*bt(2,3)+xpp(3)*bt(3,3))

       xpp(:)=0.0

       xpp(:) = xc(:)

       return
  end subroutine crist_to_cart

