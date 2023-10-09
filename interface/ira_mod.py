
# Copyright (C) 2023, MAMMASMIAS Consortium
# Written by: Miha Gunde
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: Apache-2.0
# See the file LICENSE.txt for further information.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Python interface to the library_ira.f90 and library_sofi.f90 routines, via ctypes module

from ctypes import *
import numpy as np
from os.path import dirname,abspath,join
from inspect import getsourcefile

class algo():
    def __init__(self, shlib='' ):
        # path to this file
        mypath=dirname(abspath(getsourcefile(lambda:0)))
        # one directory up
        mypath=dirname(mypath)
        # by default, the lib should be there:
        path = join(mypath,"src/libira.so")

        if shlib:
            # user provde path
            path=shlib

        self.lib = CDLL(path)

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
    designed for matching generic atomic structures.
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
            [c_int, POINTER(c_int), POINTER(c_double), \
             c_int, POINTER(c_int), POINTER(c_double), \
             c_double, \
             POINTER(POINTER(c_int*nat2)), POINTER(POINTER(c_double*nat2))]
        self.lib.lib_cshda.restype=None
        # pbc
        self.lib.lib_cshda_pbc.argtypes = \
            [c_int, POINTER(c_int), POINTER(c_double), \
             c_int, POINTER(c_int), POINTER(c_double), \
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

        The Iterative Rotatitions and Assignments (IRA) procedure to match two structures, including
        the SVD correction at the end.

        This is a wrapper to the lib_match routine in library_ira.f90

        Solution is to be applied to struc2, as:

            >>> for i in range(nat2):
            >>>     coords2[i] = np.matmul( rotation, coords2[i] ) + translation

        or alternatively to struc1 as:

            >>> for i in range(nat1):
            >>>     coords1[i] = np.matmul( rotation.T, coords1[i] ) - np.matmul( rotation.T, translation )

        Apply permutation:

            >>> coords2[:] = coords2[permutation]
            >>> typ2[:] = typ2[permutation]

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
                            should be > 1.0, for very non-congruent structures a higher value
                            is needed. Lower value speeds up the algorithm, but can give
                            non-optimal results.
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

        # manage the candidates: array of integers with proper size are needed.
        # Also shift everything by +1 to get F indices.
        # In case of equal number of atoms, default behaviour is to specify
        # signal -1 (which by default means: use geometric center); in case
        # of non-equal number of atoms, default is to take first atom as candidate
        # in structure1, and all atoms of same typ as candidate1 as candidates for structure2
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
                idx_1=candidate1+1
            else:
                # candidate1 is an array
                idx_1=candidate1[0]+1
                for i in range(lenc1):
                    cand1[i] = candidate1[i]+1
        if candidate2 is None:
            # candiate is not set,
            # matching strucs with different nat, impose all atoms of struc2
            if nat2 != nat1:
                m=0
                for i in range(nat2):
                    # skip atoms of dirrerent typ,
                    # idx_1 is atom chosen as candidate1
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

        # print(cand1)
        # print(cand2)
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

        self.lib.lib_match.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), \
              c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), \
              c_double, POINTER(POINTER(c_double*9)), POINTER(POINTER(c_double*3)), \
              POINTER(POINTER(c_int*nat2)), POINTER(c_double) ]
        self.lib.lib_match.restype=None

        self.lib.lib_match( n1, t1, c1, cd1, n2, t2, c2, cd2, \
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


class sym_data():
    """
    Definition of the object, which is returned by `ira_mod.SOFI.compute()`

    Members of the object:

    :param n_sym: number of symmetry elements found
    :type n_sym: int

    :param matrix: the list of matrices of symmetry operations
    :type matrix: np.array((n_sym, 3, 3), dtype = float)

    :param perm: the atomic permutation corresponding to application of that symmetry, start by index 0
    :type perm: np.array((n_sym, nat), dtype = int)

    :param op: Schoenflies notation: Op n^p; Id=identity, I=inversion, C=rotation, S=(roto-)inversion
    :type op: np.array((n_sym), type="U2"); size-2 strings

    :param n: Schoenflies notation n
    :type n: np.array(n_sym, dtype = int )

    :param p: Schoenflies notation p
    :type p: np.array(n_sym, dtype = int )

    :param axis: the axis over which a symmetry element operates (for S: normal vector to plane)
    :type axis: np.array((n_sym, 3), dtype = float )

    :param angle: angle of rotation of a symmetry element, in units of 1/2pi, e.g. angle=0.25 is 1/4 circle
    :type angle: np.array(n_sym, dtype = float )

    :param dmax: maximal displacement of any atom upon transformation by that symmetry matrix
    :type dmax: np.array(n_sym, dtype = float )

    :param pg: point group of the structure
    :type pg: str

    :function print: prints the full data of the `sym_data` object

    """
    def __init__(self):
        # define object to hold the resulting data
        self.n_sym=0
        self.matrix=None
        self.perm=None
        self.op=None
        self.n=None
        self.p=None
        self.axis=None
        self.angle=None
        self.dmax=None
        self.pg=None

    def print(self):
        """
        print the members of sym_data() instance
        """
        print( "sym object contains: %i symmetry operations, listed below:" % (self.n_sym) )
        for i in range(self.n_sym):
            print(i)
            print( "operation: %s %i ^ %i" % (self.op[i], self.n[i], self.p[i]) )
            print( "angle: %.4f" % self.angle[i] )
            print( "axis: %8.4f %8.4f %8.4f" %( self.axis[i][0], self.axis[i][1], self.axis[i][2]) )
            print( "matrix:" )
            print( "%8.4f %8.4f %8.4f" % ( self.matrix[i][0][0], self.matrix[i][0][1], self.matrix[i][0][2] ) )
            print( "%8.4f %8.4f %8.4f" % ( self.matrix[i][1][0], self.matrix[i][1][1], self.matrix[i][1][2] ) )
            print( "%8.4f %8.4f %8.4f" % ( self.matrix[i][2][0], self.matrix[i][2][1], self.matrix[i][2][2] ) )
            print( "dmax  %9.4f" % self.dmax[i] )
            print( "atomic permutation:")
            print( self.perm[i] )
            print()
        print("The corresponding PG is: %s" % (self.pg) )


class SOFI(algo):
    """
    This is the python interface to the SOFI shared library file.
    It should be initialized with the path to the `libsofi.so` file:

    Example:

    >>> import ira_mod
    >>> sofi = ira_mod.SOFI( '/path/to/sofi/libsofi.so' )



    """

    # nmax is imposed by sofi_tools.f90; it gives the maximal number of symmetries
    # if you need more than than then change in sofi_tools.f90 and here, and recompile sofi
    _nmax=200

    def compute(self, nat, typ_in, coords, sym_thr ):
        """
        this is a wrapper to lib_compute_all() from library.f90
        Description: This routine computes all the relevant PG symmetry data of the input structure.

        input:
        ======

        :param nat: number of atoms in the structure
        :type nat: int

        :param typ_in: atomic types of the atoms
        :type typ_in: string or int

        :param coords: the atomic positions of the structure
        :tpye coords: np.array((nat,3), dtype=float)

        :param sym_thr: symmetry threshold value
        :type sym_thr: float

        output:
        =======

        instance of `ira_mod.sym_data()` object, which cotains the full symmetry data of the input structure.

        object `sym_data` contains the full symmetry data of the input structure.
        To access info about symmetry element *i*, get the data as: `sym_data.var[i]`, where `var`
        is one of the variables in the `sym_data` object:

            `n_sym`, `matrix`, `perm`, `op`, `n`, `p`, `axis`, `angle`, `dmax`, `pg`

        see help( ira_mod.sym_data ) for full description of these variables.


        """

        ## initialize instance of sym_data
        sym = sym_data()

        ## nmax value
        nmax=self._nmax

        # check if input types are already c_int or not
        typ = self.tf_int( typ_in )

        # input data
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        thr = c_double( sym_thr )


        # output data
        nmat = c_int()
        mat_data = (c_double*9*nmax)()
        perm_data = (c_int*nat*nmax)()
        op_data = (c_char*2*nmax)()
        n_data = (c_int*nmax)()
        p_data = (c_int*nmax)()
        ax_data = (c_double*3*nmax)()
        angle_data = (c_double*nmax)()
        dmax_data = (c_double*nmax)()
        pg = (c_char*11)()


        # have to set argtypes in here, since nat is not know in init
        self.lib.lib_compute_all.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), c_double, \
              POINTER(c_int), POINTER(POINTER(c_double*9*nmax)), POINTER(POINTER(c_int*nat*nmax)), \
              POINTER(POINTER(c_char*2*nmax)), POINTER(POINTER(c_int*nmax)), POINTER(POINTER(c_int*nmax)), \
              POINTER(POINTER(c_double*3*nmax)), POINTER(POINTER(c_double*nmax)), POINTER(POINTER(c_double*nmax)), \
              POINTER(POINTER(c_char*11)) ]
        self.lib.lib_compute_all.restype=None

        # call routine from library.f90
        self.lib.lib_compute_all( n, t, c, thr, \
                                  nmat, pointer(mat_data), pointer(perm_data), \
                                  pointer(op_data), pointer(n_data), pointer(p_data), \
                                  pointer(ax_data), pointer(angle_data), pointer(dmax_data), pointer(pg) )

        # cast the result into readable things
        sym.pg = pg.value.decode()
        sym.n_sym = nmat.value
        n_op=sym.n_sym
        # set matrices, attention: order is C
        sym.matrix=np.zeros( (n_op,3,3), dtype=float)
        for n in range(n_op):
            m=0
            for i in range(3):
                for j in range(3):
                    sym.matrix[n][i][j] = mat_data[n][m]
                    m+=1

        sym.perm = np.ones((n_op, nat),dtype=int)
        for n in range(n_op):
            sym.perm[n]=perm_data[n][:]

        sym.op = np.empty(n_op,dtype="U2")
        for n in range(n_op):
            sym.op[n] = op_data[n][:].decode()

        sym.n=np.zeros(n_op,dtype=int)
        sym.n = n_data[:n_op]

        sym.p=np.zeros(n_op, dtype=int)
        sym.p = p_data[:n_op]

        sym.axis=np.zeros((n_op,3), dtype=float)
        for n in range(n_op):
            sym.axis[n] = ax_data[n][:]

        sym.angle=np.zeros(n_op, dtype=float)
        for n in range(n_op):
            sym.angle[n]=angle_data[n]

        sym.dmax=np.zeros(n_op, dtype=float)
        for n in range(n_op):
            sym.dmax[n]=dmax_data[n]
        # sym.dmax = dmax_data[:n_op]

        # return the instance of sym_data
        return sym


    def get_symm_ops(self, nat, typ_in, coords, sym_thr ):
        """
        Wrapper to the lib_get_symm_ops routine from library.f90.
        Description: finds the symmetry operation maytices of the structure in input,
        which fit the threshold `sym_thr`.

        input:
        ======

        :param nat: number of atoms in the structure
        :type nat: int

        :param typ_in: atomic types of the atoms
        :type typ_in: string or int

        :param coords: the atomic positions of the structure
        :tpye coords: np.array((nat,3), dtype=float)

        :param sym_thr: symmetry threshold value
        :type sym_thr: float


        output:
        =======

        :param n_sym: number of symmetry elements found
        :type n_sym: int

        :param matrix: the list of matrices of symmetry operations
        :type matrix: np.array((n_sym, 3, 3), dtype = float)

        """

        nmax=self._nmax

        # check if input types are already c_int or not
        typ = self.tf_int( typ_in )

        # input data
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        thr = c_double( sym_thr )


        # output data
        nmat = c_int()
        mat_data = (c_double*9*nmax)()
        # have to set argtypes in here, since nat is not know in init
        self.lib.lib_get_symm_ops.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), c_double, \
              POINTER(c_int), POINTER(POINTER(c_double*9*nmax)) ]
        self.lib.lib_get_symm_ops.restype=None

        self.lib.lib_get_symm_ops( n, t, c, thr, nmat, pointer(mat_data) )

        # cast the result into readable things
        n_op = nmat.value
        # set matrices
        matrix=np.ndarray( (n_op,3,3), dtype=float, order='C')
        for n in range(n_op):
            m=0
            for i in range(3):
                for j in range(3):
                    matrix[n][i][j] = mat_data[n][m]
                    m+=1
        return n_op, matrix


    def get_pg( self, nm_in, mat_list, verb=False):

        """
        wrapper to lib_get_pg() from library.f90
        Description: find the PG of input list of operations

        input:
        ======

        :param nm_in: number of matrices in the `mat_list` array
        :type nm_in: integer

        :param mat_list: list of matrices containging symmetry operations
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )

        == optional ==
        :param verb: verbosity of the output, verb=True outputs a report from get_pg() routine.
        :type verb: logical


        output:
        =======

        :param pg: associated Point Group (PG)
        :type pg: string

        """
        # input data
        n = c_int( nm_in )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )
        cverb = c_bool( verb )


        # output data
        pg = (c_char*11)()

        # have to set argtypes in here, since nat is not know in init
        self.lib.lib_get_pg.argtypes = \
            [ c_int, POINTER(c_double), POINTER(POINTER(c_char*11)), c_bool ]
        self.lib.lib_get_pg.restype=None
        self.lib.lib_get_pg( n, mats, pointer(pg), cverb)

        pg=pg.value.decode()
        return pg

    def get_unique_ax_angle( self, nm_in, mat_list ):
        """
        wrapper to lib_get_unique_ax_angle() from library.f90
        Description: input list of symmetry matrices, output lists
        of op, ax, angle for each symmetry operation in the list, such that
        axes are aligned, and angle has a correspoding +/- sign.

        This is a wrapper to lib_unique_ax_angle() routine from library_sofi.f90

        input:
        ======

        :param nm_in: number of matrices in the `mat_list` array
        :type nm_in: integer

        :param mat_list: list of matrices containging symmetry operations
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )


        output:
        =======

        :param op: Schoeflies Op name, possible values:
                         - "Id" = identity
                         - "I" = inversion
                         - "C" = rotation
                         - "S" = (roto-)reflection
        :type op: array of length-2 strings, nd.ndarray( nm_in, dtype = "U2" )

        :param axis: the list of axes along each corresponding matrix acts.
        :type axis: nm_in x 3 vector, np.ndarray((nm_in, 3), dtype = float)

        :param angle: the list of angles by which each matrix rotates, in units of 1/(2pi),
                      e.g. the value `angle = 0.5` means half the circle, and `angle = -1.0/6.0`
                      means rotation about the axis in negative direction for 1/6th of the circle
        :type angle: np.ndarray( nm_in, dtype = float )

        """
        # input data
        n_m = c_int( nm_in )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )
        n_op=nm_in
        nl = 2*n_op+1

        # output
        op_data = (c_char*nl)()
        ax_data = (c_double*3*n_op)()
        angle_data = (c_double*n_op)()

        self.lib.lib_unique_ax_angle.argtypes = \
            [ c_int, POINTER(c_double), POINTER(POINTER(c_char*nl)) , \
              POINTER(POINTER(c_double*3*n_op)), POINTER(POINTER(c_double*n_op)) ]
        self.lib.lib_unique_ax_angle.restype=None

        self.lib.lib_unique_ax_angle( n_m, mats, pointer(op_data), \
                                      pointer(ax_data), pointer(angle_data) )
        op = np.empty(n_op,dtype="U2")
        for n in range(n_op):
            op[n] = op_data[2*n:2*n+2].decode()
        axis=np.ndarray((n_op,3), dtype=float)
        for n in range(n_op):
            axis[n] = ax_data[n][:]
        angle=np.zeros(n_op, dtype=float)
        for n in range(n_op):
            angle[n]=angle_data[n]
        # angle=angle_data

        return op, axis, angle


    def analmat( self, rmat ):
        """
        wrapper to lib_analmat.f90 from library.f90
        Description: analyse the input matrix, output the Schoenflies notation: Op n^p
        The notation can be transformed into an angle:

            angle = p / n * 2pi

        e.g. C 8^3 corresponds to 3/8th of full circle = 135 degress = 2.3563 radian

        Also return axis and angle of the matrix

        This is a wrapper to lib_analmat() routine from library_sofi.f90

        input:
        ======

        :param rmat: input matrix
        :type rmat: np.ndarray((3, 3), dtype = float )


        output:
        =======

        :param op: Schoeflies Op name, possible values:
                         - "Id" = identity
                         - "I" = inversion
                         - "C" = rotation
                         - "S" = (roto-)reflection
        :type op: string

        :param n: the n value from Op n^p, gives the angle 2pi/n
        :type n: integer

        :param p: the p value from Op n^p, gives the "multiplicity" of 1/n operation
        :type p: integer

        :param ax: axis along which the matrix `rmat` operates, in case of rotation is the axis of rotation,
                   in case of (roto-)reflection is the normal of the plane of reflection.
        :type ax: 3D vector, np.ndarray(3, dtype = float)

        :param angle: the angle of rotation of matrix `rmat`, in units of 1/2pi, is the ratio p/n
        :type angle: float

        """

        # input data
        mat = rmat.ctypes.data_as( POINTER(c_double) )

        # output
        op=(c_char*3)()
        n=c_int()
        p=c_int()
        ax=(c_double*3)()
        angle=c_double()

        self.lib.lib_analmat.argtypes=\
            [ POINTER(c_double), POINTER(POINTER(c_char*3)), \
              POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_double*3)), \
              POINTER(c_double) ]
        self.lib.lib_analmat.restype=None

        self.lib.lib_analmat( mat, pointer(op), pointer(n), pointer(p), \
                              pointer(ax), pointer(angle) )


        op=op.value.decode()
        n=n.value
        p=p.value
        angle=angle.value
        ax=ax[:]

        return op, n, p, ax, angle


    def ext_Bfield( self, n_mat, mat_list, b_field ):

        """
        wrapper to lib_ext_Bfield() from library.f90
        Description: Impose a constraint on the relevant symmetry operations,
        in the form of an external magnetic field as arbitrary vector.

        Implementation of the functionality proposed in:


            Ansgar Pausch, Melanie Gebele, Wim Klopper;
            Molecular point groups and symmetry in external magnetic fields.
            J. Chem. Phys. 28 November 2021; 155 (20): 201101.
            DOI: https://doi.org/10.1063/5.0069859


        The magnetic field acts as an "external constraint" on the symmetry matrices,
        and can thus be seen as a post-processing of the unconstrained (full) PG of a structure.

        Input:
        ======

        :param n_mat: number of matrices in the input `mat_list`
        :type n_mat: integer

        :param mat_list: list of 3x3 matrices of symmetry operations in input
        :type mat_list: np.ndarray((n_mat, 3, 3), dtype = float )

        :param b_field: direction of magnetic field
        :type b_field: 3d vector, np.ndarray(3, dtype = float )


        output:
        =======

        :param n_op: number of operations that fulfill the constraint of the desired magnetic field
        :type n_op: integer

        :param matrix: the list of matrices that act on the system constrained by the desired magnetic field
        :type matrix: np.ndarray( (n_op, 3, 3), dtype = float )

        """
        # input data
        n_m = c_int( n_mat )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )
        bf = b_field.ctypes.data_as( POINTER(c_double) )
        # output
        n_out=c_int()
        m_out=(c_double*9*n_mat)()

        self.lib.lib_ext_bfield.argtypes = [ c_int, POINTER(c_double), \
                                             POINTER(c_double), POINTER(c_int), \
                                             POINTER(POINTER(c_double*9*n_mat))]
        self.lib.lib_ext_bfield.restype=None

        self.lib.lib_ext_bfield( n_m, pointer(mats), bf, pointer(n_out), pointer(m_out) )

        n_op=n_out.value
        matrix=np.ndarray( (n_op,3,3), dtype=float)
        for n in range(n_op):
            m=0
            for i in range(3):
                for j in range(3):
                    matrix[n][i][j] = m_out[n][m]
                    m+=1

        return n_op, matrix


    def get_perm( self, nat, typ_in, coords, nm_in, mat_list ):

        '''

        Description:
        Obtain the list of permutations and maximal distances for each matrix in input.

        This is a wrapper to lib_get_perm() from library_sofi.f90

        Input:
        ======

        :param nat: number of atoms in the structure
        :type nat: integer

        :param typ_in: atomic types of the structure
        :type typ_in: np.ndarray( nat, dtype = integer or string )

        :param coords: atomic positions in the structure
        :type coords: np.ndarray((nat, 3), dtype = float )

        :param nm_in: number of matrices in the `mat_list` input
        :typ nm_in: integer

        :param mat_list: list of input matrices
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )

        Output:
        =======

        :param perm: list of permutations corresponding to each matrix operation in input
        :type perm: np.ndarray( (nm_in, nat), dtype = integer )

        :param dmax: maximal distance (Hausdorff distance) corresponding to the application off
                     each matrix in the `mat_list` to the structure. Can be seen as "score" of each symmetry
        :type dmax: np.ndarray( nm_in, dtype = float )

        '''

        # check if input types are already c_int or not
        typ = self.tf_int( typ_in )

        # input data
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        n_m = c_int( nm_in )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )

        # output data
        perm_data = (c_int*nat*nm_in)()
        dmax_data = (c_double*nm_in)()


        # have to set argtypes in here, since nat is not know in init
        self.lib.lib_get_perm.argtypes = \
            [ c_int, \
              POINTER(c_int), \
              POINTER(c_double), \
              c_int, \
              POINTER(c_double), \
              POINTER(POINTER(c_int*nat*nm_in)), \
              POINTER(POINTER(c_double*nm_in)) ]
        self.lib.lib_get_perm.restype=None

        self.lib.lib_get_perm( n, t, c, n_m, mats, \
                               pointer(perm_data), pointer(dmax_data) )

        # cast the result into readable things
        perm = np.ndarray((nm_in, nat),dtype=int)
        for n in range(nm_in):
            perm[n]=perm_data[n][:]
        dmax=np.ndarray(nm_in, dtype=float)
        for n in range(nm_in):
            dmax[n] = dmax_data[n]

        return perm, dmax

    def get_combos( self, nat, typ_in, coords, nm_in, mat_list ):

        '''
        Description:
        Obtain all unique combinations of input matrices, which are symmetry operations of given structure.

        This is a wrapper to the routine lib_get_combos() from library_sofi.f90

        Input:
        ======

        :param nat: number of atoms in the structure
        :type nat: integer

        :param typ_in: atomic types of the structure
        :type typ_in: np.ndarray( nat, dtype = integer or string )

        :param coords: atomic positions in the structure
        :type coords: np.ndarray((nat, 3), dtype = float )

        :param nm_in: number of matrices in the `mat_list` input
        :typ nm_in: integer

        :param mat_list: list of input matrices
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )


        Output:
        =======
        :param nm_out: number of matrices in the output list `mat_out`
        :type nm_out: integer

        :param mat_out: list of output matrices
        :type mat_out: np.ndarray( (nm_out, 3, 3), dtype = float )

        '''

        nmax=self._nmax
        # check if input types are already c_int or not
        typ = self.tf_int( typ_in )

        #input
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        nm = c_int( nm_in )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )

        # output
        nb_out= c_int()
        cmat_out=(c_double*9*nmax)()

        self.lib.lib_get_combos.argtypes=\
            [ c_int, \
              POINTER(c_int), \
              POINTER(c_double), \
              c_int, \
              POINTER(c_double), \
              POINTER(c_int), \
              POINTER(POINTER(c_double*9*nmax)) ]
        self.lib.lib_get_combos.restype=None

        self.lib.lib_get_combos( n, t, c, nm, mats, nb_out, pointer(cmat_out) )

        nm_out=nb_out.value
        mat_out=np.ndarray( (nm_out,3,3), dtype=float)
        for n in range(nm_out):
            m=0
            for i in range(3):
                for j in range(3):
                    mat_out[n][i][j] = cmat_out[n][m]
                    m+=1
        return nm_out, mat_out


    def try_mat( self, nat, typ_in, coords, theta ):
        '''
        Description:
        Apply a given matrix `theta` to the structure, and compute the distance between the transformed
        and original structure. If the distance is small, then `theta` is a symmetry operation (or close to).
        The distance is computed in a permutationally invariant way, using the CShDA algorithm, which
        imposes a one-to-one assignment of the atoms. The value of distance corresponds to the maximal value
        of distances between all assigned atomic pairs, which is the Hausdorff distance.

        This is a wrapper to lib_try_mat() routine from library_sofi.f90

        Input:
        ======

        :param nat: number of atoms in the structure
        :type nat: integer

        :param typ_in: atomic types of the structure
        :type typ_in: np.ndarray( nat, dtype = integer or string )

        :param coords: atomic positions in the structure
        :type coords: np.ndarray((nat, 3), dtype = float )

        :param theta: the input matrix of an orthonormal operation to be tested on structure
        :type theta: np.ndarray((3,3), dtype = float)


        Output:
        =======

        :param dh: the Hausdorff distance value
        :type dh: float

        :param perm: permutations of atoms which were assignmed by CShDA upon application of `theta`
        :type perm: np.ndarray((nat), dtype = int )

        '''
        typ = self.tf_int( typ_in )

        #input
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        r = theta.ctypes.data_as( POINTER(c_double) )

        # output
        c_dh=c_double()
        c_perm=(c_int*nat)()

        self.lib.lib_try_mat.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), \
              POINTER(c_double), POINTER(c_double), POINTER(POINTER(c_int*nat)) ]
        self.lib.lib_try_mat.restype = None

        self.lib.lib_try_mat( n, t, c, r, pointer(c_dh), pointer(c_perm) )


        dh = c_dh.value
        perm = np.ndarray(nat, dtype=int)
        for i in range(nat):
            perm[i]=c_perm[i]

        return dh, perm

    def construct_operation( self, op, axis, angle ):
        '''
        Description:
        Construct the 3x3 matrix corresponding to operation encoded by the Schoenflies Op,
        such that it acts along the desired axis, and given angle.

        This is a wrapper to lib_construct_operation() routine from library_sofi.f90

        Input:
        ======

        :param op: Schoeflies Op name, possible values:
                         - "Id" = identity
                         - "I" = inversion
                         - "C" = rotation
                         - "S" = (roto-)reflection
        :type op: string

        :param axis: the axis along which the desired operation should act
        :type axis: 3D vector, np.ndarray(3, dtype = float)

        :param angle: the angle by which the desired operation should rotate, in units of 1/(2pi),
                      e.g. the value `angle = 0.5` means half the circle, and `angle = -1.0/6.0`
                      means rotation about the axis in negative direction for 1/6th of the circle
        :type angle: float


        Output:
        =======

        :param matrix: the matrix representation of desired operation
        :type matrix: 3x3 matrix, np.ndarray((3, 3), dtype = float )


        NOTE:
        =====

        The operation "Id" is equivalent to "C" about any axis for angle=0.0; and
        similarly the operation "I" is equivalent to "S" with angle=0.5 about any axis.
        The operation "S" with angle=0.0 is a reflection.

        '''
        c_op=op.encode()
        c_ax = axis.ctypes.data_as(POINTER(c_double))
        c_angle = c_double(angle)

        mat_out=(c_double*9)()

        self.lib.lib_construct_operation.argtypes = \
            [ c_char_p, POINTER(c_double), c_double, POINTER(POINTER(c_double*9)) ]
        self.lib.lib_construct_operation.restype = None

        self.lib.lib_construct_operation( c_op, c_ax, c_angle, pointer(mat_out) )

        matrix = np.ndarray((3,3), dtype=float)
        m=0
        for i in range(3):
            for j in range(3):
                matrix[i][j] = mat_out[m]
                m+=1
        return matrix

    def mat_combos( self, nm_in, mat_list ):

        '''
        Description:
        Obtain all unique combinations of matrices in input, without checking them against any
        specific structure.

        This is a wrapper to routine lib_mat_combos() from library_sofi.f90

        Input:
        ======
        :param nm_in: number of matrices in the `mat_list` input
        :typ nm_in: integer

        :param mat_list: list of input matrices
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )


        Output:
        =======
        :param nm_out: number of matrices in the output list `mat_out`
        :type nm_out: integer

        :param mat_out: list of output matrices
        :type mat_out: np.ndarray( (nm_out, 3, 3), dtype = float )

        '''

        nmax=self._nmax

        #input
        nm = c_int( nm_in )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )

        # output
        nb_out= c_int()
        cmat_out=(c_double*9*nmax)()

        self.lib.lib_mat_combos.argtypes=\
            [ \
              c_int, \
              POINTER(c_double), \
              POINTER(c_int), \
              POINTER(POINTER(c_double*9*nmax)) \
             ]
        self.lib.lib_mat_combos.restype=None

        self.lib.lib_mat_combos( nm, mats, nb_out, pointer(cmat_out) )

        nm_out=nb_out.value
        mat_out=np.ndarray( (nm_out,3,3), dtype=float)
        for n in range(nm_out):
            m=0
            for i in range(3):
                for j in range(3):
                    mat_out[n][i][j] = cmat_out[n][m]
                    m+=1
        return nm_out, mat_out

