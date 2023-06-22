
# Copyright (C) 2023, MAMMASMIAS Consortium
# Written by: Miha Gunde

# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: Apache-2.0
# See the file LICENSE.txt for further information.

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Python interface to the library_ira.f90 routines, via ctypes module

from ctypes import *
import numpy as np


class algo():
    def __init__(self, shlib=None ):
        if shlib is None:
            raise ValueError("Please provide path to shared library file shlib_ira.so")
        else:
            self.lib = CDLL(shlib)

    def tf_int(self, arr_in):
        """
        function to transform single array arr_in into array of c_int with unqiue values.
        """
        u, indices = np.unique( arr_in, return_inverse=True )
        arr_out = np.intc( indices+1 )
        return arr_out

    def tf_int2( self, arr1_in, arr2_in ):
        """
        function to transform two integer arrays `arr1_in` and `arr2_in` into arrays of c_int
        with simlutaneously unique values
        """
        ## get list of unique values
        u=np.unique(arr1_in)
        arr1_out=np.ndarray(len(arr1_in), dtype=int)
        arr2_out=np.ndarray(len(arr2_in), dtype=int)
        ## find values according to u
        for i in range(len(arr1_in)):
            loc = np.where( arr1_in[i] == u)[0]
            if len(loc) > 0:
                arr1_out[i] = loc[0] + 1
        for i in range(len(arr2_in)):
            loc = np.where( arr2_in[i] == u)[0]
            if len(loc) > 0:
                arr2_out[i] = loc[0] + 1
            else:
                ## it can happen that arr2 has values which arr1 doesn't

                ## append the new typ into u
                u=np.append(u, arr2_in[i] )
                loc = np.where( arr2_in[i] == u)[0]
                arr2_out[i] = loc[0]+1

        arr1_out = np.intc( arr1_out )
        arr2_out = np.intc( arr2_out )

        return arr1_out, arr2_out

class IRA(algo):
    """
    The Iterative Rotatitions and Assignments (IRA) algorithm is a shape matching algorithm
    designed for matching arbitrary atomic structures.
    It solves the problem:

        P_B B = R A + t

    where `A` and `B` are two atomic structures, `R` is a rotation matrix, `t` is a translation vector,
    and `P_B` is a permutation matrix.

    For futher details please see the publication:

    M Gunde, N Salles, A Hemeryck, L Martin Samos:
    J. Chem. Inf. Model. 2021, 61, 11, 5446–5457
    DOI: https://doi.org/10.1021/acs.jcim.1c00567
    arXiv: https://export.arxiv.org/abs/2111.00939

    or the thesis:
    M Gunde: "Development of IRA : a shape matching algorithm, its implementation,
    and utility in a general off-lattice kMC kernel", 2021
    url: https://theses.hal.science/tel-03635139/

    ================================================================================


    """

    def cshda( self, nat1, typ1, coords1, nat2, typ2, coords2, lat=None, thr=None):
        """

        Constrained Shortest Distance Assignment (CShDA) algorithm, is a LAP
        (Linear Assignment Problem) solver, which constrains the atoms to be assigned
        such that distances are minimised locally. This is sometimes referred to as
        the "inverse bottleneck assignment".

        CShDA can be used to compute the assignment of atoms and the corresponding
        distances between atoms in permuted structures.
        It is used inside the IRA algorithm.

        ===========================================

        This is a wrapper to the lib_cshda and lib_cshda_pbc routines from library_ira.f90

        Requirement: nat1 .le. nat2 !!

        input:
        ======

        :param nat1: number of atoms of structure 1
        :type nat1: int

        :param typ1: atomic types of structure 1
        :type typ1: np.ndarray( nat1, dtype = int or string )

        :param coords1: atomic positions of structure 1
        :type coords1: np.ndarray( (nat1, 3), dtype = float )

        :param nat2: number of atoms of structure 2
        :type nat2: int

        :type typ2: atomic types of structure 2
        :param typ2: np.ndarray( nat2, dtype = int or string )

        :param coords2: atomic positions of structure 2
        :type coords2: np.ndarray( (nat2, 3), dtype = float)

        == optional ==

        :param lat: matrix of lattice vectors in rows. If not None, it triggers `cshda_pbc` routine.
                    Use with care.
        :type lat: np.ndarray( (3, 3), dtype = float )

        :param thr: threshold for early exit of cshda. Use with care.
        :type thr: float


        :raises ValueError: when the requirement `nat1 =< nat2` is not met.


        output:
        =======

        :param found: array of assignments of atoms of structure1 to structure2,
                      e.g.: `found[3]=4` means the atom index 3 from structure1 has been assigned to
                      atom index 4 from structure2. When `nat1 != nat2`, the first nat1 indices
                      give the permutation, the rest stay fixed.
        :type found: np.ndarray( nat2, dtype = int )

        :param dists: array of distances from atom index `i` in structure1 to atom `found[i]` in structure2.
                      The Hausdorff distance can be obtained as `dH = np.max( dists )`
        :type dists: np.ndarray( nat1, dtype = float )

        """

        if nat1 > nat2:
            error_msg='''Number of atoms nat1 must be less or equal to nat2!
            You can switch the labelling: struc1 <==> struc2.'''
            raise ValueError( error_msg )

        # check typ arrays
        u_typ1, u_typ2  = self.tf_int2( typ1, typ2 )

        # convert input to C-style
        n1 = c_int(nat1)
        n2 = c_int(nat2)
        t1 = u_typ1.ctypes.data_as( POINTER(c_int) )
        t2 = u_typ2.ctypes.data_as( POINTER(c_int) )
        c1 = coords1.ctypes.data_as( POINTER(c_double) )
        c2 = coords2.ctypes.data_as( POINTER(c_double) )
        if lat is not None:
            l = lat.ctypes.data_as( POINTER(c_double) )
        cthr = c_double(99.9)
        if thr is not None:
            cthr = c_double(thr)

        # allocate output arrays
        cfound = (c_int*nat2)()
        cdists = (c_double*nat2)()

        # non-pbc
        self.lib.lib_cshda.argtypes = \
            [c_int, POINTER(c_int), POINTER(POINTER(c_double)), \
             c_int, POINTER(c_int), POINTER(POINTER(c_double)), \
             c_double, \
             POINTER(POINTER(c_int*nat2)), POINTER(POINTER(c_double*nat2))]
        self.lib.lib_cshda.restype=None
        # pbc
        self.lib.lib_cshda_pbc.argtypes = \
            [c_int, POINTER(c_int), POINTER(POINTER(c_double)), \
             c_int, POINTER(c_int), POINTER(POINTER(c_double)), \
             POINTER(c_double), c_double,  \
             POINTER(POINTER(c_int*nat2)), POINTER(POINTER(c_double*nat2))]
        self.lib.lib_cshda_pbc.restype=None

        if lat is None:
            # call non-pbc
            self.lib.lib_cshda( n1, t1, c1, n2, t2, c2, cthr, pointer(cfound), pointer(cdists) )
        else:
            # call pbc
            self.lib.lib_cshda_pbc( n1, t1, c1, n2, t2, c2, l, cthr, pointer(cfound), pointer(cdists) )

        ## output
        found=np.ndarray(nat2, dtype=int)
        dists=np.ndarray(nat2, dtype=float)
        for i in range(nat2):
            found[i] = cfound[i]
            dists[i] = cdists[i]

        return found, dists


    def match( self, nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor, candidate1=None, candidate2=None ):
        """

        The Iterative Rotatitions and Assignments (IRA) procedure to match two structures.

        This is a wrapper to the lib_ira_unify routine in library_ira.f90

        Solution is to be applied to struc2, as:

            >>> for i in range(nat2):
            >>>     coords2[i] = np.matmul( rotation, coords2[i] ) + translation

        or alternatively to struc1 as:

            >>> for i in range(nat1):
            >>>     coords1[i] = np.matmul( rotation.T, coords1[i] ) - np.matmul( rotation.T, translation )

        Apply permutation:

            >>> coords2[:] = coords2[permutation]

        Requirement: nat1 .le. nat2 !!
        ==========================

        input:
        ======
        :param nat1: number of atoms of structure 1
        :type nat1: int

        :param typ1: atomic types of structure 1
        :type typ1: np.ndarray( nat1, dtype = int or string )

        :param coords1: atomic positions of structure 1
        :type coords1: np.ndarray( (nat1, 3), dtype = float )

        :param nat2: number of atoms of structure 2
        :type nat2: int

        :type typ2: atomic types of structure 2
        :param typ2: np.ndarray( nat2, dtype = int or string )

        :param coords2: atomic positions of structure 2
        :type coords2: np.ndarray( (nat2, 3), dtype = float)

        :param kmax_factor: factor for multiplication of search radius, the value
                            should be > 1.0, for very non-congruent structures a higher value is needed.
                            Lower value speeds up the algorithm, but can give non-optimal values.
        :type kmax_factor: float

        === optional ===
        :param candidate1: list of candidate central atoms of structure1. This is useful when some
                           additional information is known about the structures, e.g. we want to impose
                           some specific atoms as possible central atoms. Particularly useful when
                           matching structures with different number of atoms.
                           If None, geometric center is used as central point.
                           NOTE: the indices should be starting from 0. Use with care.
        :type candidate1: int, or np.ndarray( AnySize, dtype = int )

        :param candidate2: similar as candidate1, except for structure2
        :type candidate2: int, or np.ndarray( AnySize, dtype = int )

        :raises ValueError: when kmax_factor is =< 1.0
        :raises ValueError: when the requirement `nat1 =< nat2` is not met.


        output:
        =======

        :param rotation: rotation matrix
        :type rotation: np.ndarray( (3,3), dtype = float)

        :param translation: translation vector
        :type translation: np.ndarray( 3, dtype = float)

        :param permutation: permutation of atoms
        :type permutation: np.ndarray( (nat2), dtype = int )

        :param hd: Hausdorff distance value of the match
        :type hd: float

        """

        if nat1 > nat2:
            error_msg='''Number of atoms nat1 must be less or equal to nat2!
            You can switch the labelling: struc1 <==> struc2.'''
            raise ValueError( error_msg )

        if kmax_factor <= 1.0:
            error_msg='''The kmax_factor should have a value larger than 1.0!
            Larger value is needed for non-congruent structures.'''
            raise ValueError( error_msg )

        # check typ arrays
        u_typ1, u_typ2  = self.tf_int2( typ1, typ2 )

        # manage the candidates: array of integers with proper size are needed
        # also shift everything by +1 to get F indices
        cand1 = np.zeros(nat1, dtype=int)
        cand2 = np.zeros(nat2, dtype=int)
        idx_1=0
        if candidate1 is None:
            # candidate1 is not set
            if nat2 != nat1:
                # matching strucs with different nat, impose first atom of struc1
                cand1[0]=idx_1+1
            else:
                # matching equal nat, signal -1
                cand1[0]=-1
        else:
            lenc1 = np.size(candidate1)
            if lenc1 == 1:
                # candidate1 is a single index
                cand1[0]=candidate1+1
            else:
                # candidate1 is an array
                for i in range(lenc1):
                    cand1[i] = candidate1[i]+1
        if candidate2 is None:
            # candiate is not set

            # matching strucs with different nat, impose all atoms of struc2
            if nat2 != nat1:
                m=0
                for i in range(nat2):
                    if typ2[i] != typ1[idx_1]:
                        continue
                    cand2[m]=i+1
                    m+=1
            # matching equal nat, signal -1
            else:
                cand2[0]=-1
        else:
            lenc2 = np.size(candidate2)
            if lenc2 == 1:
                # candidate2 is a single index
                cand2[0]=candidate2+1
            else:
                # candidate2 is an array
                for i in range(lenc2):
                    cand2[i] = candidate2[i]+1

        # convert input to C-style
        n1 = c_int(nat1)
        n2 = c_int(nat2)
        t1 = u_typ1.ctypes.data_as( POINTER(c_int) )
        t2 = u_typ2.ctypes.data_as( POINTER(c_int) )
        c1 = coords1.ctypes.data_as( POINTER(c_double) )
        c2 = coords2.ctypes.data_as( POINTER(c_double) )
        cd1 = np.intc(cand1).ctypes.data_as( POINTER(c_int) )
        cd2 = np.intc(cand2).ctypes.data_as( POINTER(c_int) )
        km = c_double( kmax_factor )

        # allocate C output
        c_rmat = (c_double*9)()
        c_tr = (c_double*3)()
        c_perm = (c_int*nat2)()
        c_hd = c_double()

        self.lib.lib_ira_unify.argtypes = \
            [ c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_int), \
              c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_int), \
              c_double, POINTER(POINTER(c_double*9)), POINTER(POINTER(c_double*3)), \
              POINTER(POINTER(c_int*nat2)), POINTER(c_double) ]
        self.lib.lib_ira_unify.restype=None

        self.lib.lib_ira_unify( n1, t1, c1, cd1, n2, t2, c2, cd2, \
                                km, pointer(c_rmat), pointer(c_tr), pointer(c_perm), pointer(c_hd) )

        # convert output C data
        rotation = np.ndarray((3,3),dtype=float)
        m=0
        for i in range(3):
            for j in range(3):
                rotation[i][j] = c_rmat[m]
                m+=1

        # output
        permutation = np.ndarray(nat2, dtype=int)
        for i in range(nat2):
            permutation[i] = c_perm[i]
        translation = np.ndarray(3,dtype=float)
        for i in range(3):
            translation[i] = c_tr[i]
        hd=c_hd.value

        return rotation, translation, permutation, hd


