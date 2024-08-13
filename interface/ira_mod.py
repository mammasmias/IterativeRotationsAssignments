
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
    """
    This is the python interface to the IRA and SOFI shared library file `libira.so`.
    In order to use it, the path to this file should be added to
    the environment variable PYTHONPATH:

    .. code-block:: bash

        export PYTHONPATH=$PYTHONPATH:/your/path/to/IRA/interface


    Then IRA and SOFI can be imported and used in python as:

    >>> import ira_mod
    >>> ira = ira_mod.IRA()
    >>> sofi = ira_mod.SOFI( )



    """

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
        self.__version__, self.__date__ = self.get_version()

    def tf_int(self, arr_in):
        """
        function to transform single array arr_in into array of c_int with unqiue values.

        :meta private:
        """
        u, indices = np.unique( arr_in, return_inverse=True )
        arr_out = np.intc( indices+1 )
        return arr_out

    def tf_int2( self, arr1_in, arr2_in ):
        """
        function to transform two integer arrays `arr1_in` and `arr2_in` into arrays of c_int
        with simlutaneously unique values

        :meta private:
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

    def get_version( self ):
        '''
        obtain the IRA library version
        '''
        self.lib.libira_get_version.restype = None
        self.lib.libira_get_version.argtypes = [ c_char_p, POINTER(c_int) ]
        cstring = (c_char*10)()
        cdate = c_int()

        self.lib.libira_get_version( cstring, cdate )

        string = cstring.value.decode()
        date = cdate.value

        return string, date

class IRA(algo):
    """
    The Iterative Rotatitions and Assignments (IRA) algorithm is a shape matching algorithm
    designed for matching generic atomic structures.
    It solves the problem:

    .. math::
        P_B B = \mathbf{R} A + \mathbf{t}


    where :math:`A` and :math:`B` are two atomic structures, :math:`\mathbf{R}` is a rotation matrix,
    :math:`\mathbf{t}` is a translation vector,
    and :math:`P_B` is a permutation matrix.

    For futher details please see the publication:

    M Gunde, N Salles, A Hemeryck, L Martin Samos:
    J. Chem. Inf. Model. 2021, 61, 11, 5446–5457
    DOI: https://doi.org/10.1021/acs.jcim.1c00567
    arXiv: https://export.arxiv.org/abs/2111.00939

    or the thesis:
    M Gunde: "Development of IRA : a shape matching algorithm, its implementation,
    and utility in a general off-lattice kMC kernel", 2021
    url: https://theses.hal.science/tel-03635139/


    In order to use the module, the path to this file should be added to the environment variable PYTHONPATH:

    .. code-block:: bash

        export PYTHONPATH=/your/path/to/IRA/interface


    Example import:

    >>> import ira_mod
    >>> sofi = ira_mod.IRA( )


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

        This is a wrapper to the libira_cshda and libira_cshda_pbc routines from library_ira.f90

        .. note::
           Requirement: nat1 :math:`\le` nat2

        **== input: ==**

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

        **== optional ==**

        :param lat: matrix of lattice vectors in rows. If not None, it triggers `cshda_pbc` routine.
                    Use with care.
        :type lat: np.ndarray( (3, 3), dtype = float )

        :param thr: threshold for early exit of cshda. Use with care.
        :type thr: float


        **== output: ==**

        :param found: array of assignments of atoms of structure1 to structure2,
                      e.g.: `found[3]=4` means the atom index 3 from structure1 has been assigned to
                      atom index 4 from structure2. When `nat1 != nat2`, the first nat1 indices
                      give the permutation, the rest stay fixed.
        :type found: np.ndarray( nat2, dtype = int )

        :param dists: array of distances from atom index `i` in structure1 to atom `found[i]` in structure2.
                      The Hausdorff distance can be obtained as `dH = np.max( dists )`
        :type dists: np.ndarray( nat1, dtype = float )


        :raises ValueError: when the requirement `nat1 =< nat2` is not met.


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
        self.lib.libira_cshda.argtypes = \
            [c_int, POINTER(c_int), POINTER(c_double), \
             c_int, POINTER(c_int), POINTER(c_double), \
             c_double, \
             POINTER(POINTER(c_int*nat2)), POINTER(POINTER(c_double*nat2))]
        self.lib.libira_cshda.restype=None
        # pbc
        self.lib.libira_cshda_pbc.argtypes = \
            [c_int, POINTER(c_int), POINTER(c_double), \
             c_int, POINTER(c_int), POINTER(c_double), \
             POINTER(c_double), c_double,  \
             POINTER(POINTER(c_int*nat2)), POINTER(POINTER(c_double*nat2))]
        self.lib.libira_cshda_pbc.restype=None

        if lat is None:
            # call non-pbc
            self.lib.libira_cshda( n1, t1, c1, n2, t2, c2, cthr, pointer(cfound), pointer(cdists) )
        else:
            # call pbc
            self.lib.libira_cshda_pbc( n1, t1, c1, n2, t2, c2, l, cthr, pointer(cfound), pointer(cdists) )

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

        This function is a wrapper to the libira_match routine in library_ira.f90

        It returns the solution to:

        .. math::
            P_B B = \mathbf{R} A + \mathbf{t}


        where :math:`A` and :math:`B` are two atomic structures, :math:`\mathbf{R}` is a rotation matrix,
        :math:`\mathbf{t}` is a translation vector,
        and :math:`P_B` is a permutation matrix.

        Solution is to be applied to struc2, as:

            >>> for i in range(nat2):
            >>>     coords2[i] = np.matmul( rotation, coords2[i] ) + translation

        or alternatively to struc1 as:

            >>> for i in range(nat1):
            >>>     coords1[i] = np.matmul( rotation.T, coords1[i] ) - np.matmul( rotation.T, translation )

        Apply permutation:

            >>> coords2[:] = coords2[permutation]
            >>> typ2[:] = typ2[permutation]

        .. note::
            Requirement: nat1 :math:`\le` nat2

        **== input: ==**

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

        **== optional ==**

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


        **== output: ==**

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
        c_err = c_int()

        self.lib.libira_match.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), \
              c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), \
              c_double, POINTER(POINTER(c_double*9)), POINTER(POINTER(c_double*3)), \
              POINTER(POINTER(c_int*nat2)), POINTER(c_double), POINTER(c_int) ]
        self.lib.libira_match.restype=None

        self.lib.libira_match( n1, t1, c1, cd1, n2, t2, c2, cd2, \
                            km, pointer(c_rmat), pointer(c_tr), pointer(c_perm), pointer(c_hd), pointer(c_err) )
        if c_err.value != 0:
            msg = "error"
            raise ValueError(msg)

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
    .. _sym_data:

    Definition of the object, which is returned by :ref:`ira_mod.SOFI.compute() <sofi.compute>`

    Members of the object:

    :param n_sym: number of symmetry elements found
    :type n_sym: int

    :param matrix: the list of matrices of symmetry operations
    :type matrix: np.ndarray((n_sym, 3, 3), dtype = float)

    :param perm: the atomic permutation corresponding to application of that symmetry, start by index 0
    :type perm: np.ndarray((n_sym, nat), dtype = int)

    :param op: Schoenflies notation: Op n^p; E=identity, I=inversion, C=rotation, S=(roto-)inversion
    :type op: np.ndarray((n_sym), type="U1"); size-1 strings

    :param n: Schoenflies notation n
    :type n: np.ndarray(n_sym, dtype = int )

    :param p: Schoenflies notation p
    :type p: np.ndarray(n_sym, dtype = int )

    :param axis: the axis over which a symmetry element operates (for S: normal vector to plane)
    :type axis: np.ndarray((n_sym, 3), dtype = float )

    :param angle: angle of rotation of a symmetry element, in units of 1/2pi, e.g. angle=0.25 is 1/4 circle
    :type angle: np.ndarray(n_sym, dtype = float )

    :param dHausdorff: maximal displacement of any atom upon transformation by that symmetry matrix
    :type dHausdorff: np.ndarray(n_sym, dtype = float )

    :param pg: point group of the structure
    :type pg: str

    :param prin_ax: principal axis of the PG
    :type prin_ax: np.ndarray( 3, dtype = float )

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
        self.dHausdorff=None
        self.pg=None
        self.prin_ax=None

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
            print( "dHausdorff  %9.4f" % self.dHausdorff[i] )
            print( "atomic permutation:")
            print( self.perm[i] )
            print()
        print("The corresponding PG is: %s" % (self.pg) )
        print("List of principal axes, N = %i :" %(self.n_prin_ax) )
        for i in range( self.n_prin_ax ):
           print("%8.4f %8.4f %8.4f" % (self.prin_ax[i][0], self.prin_ax[i][1], self.prin_ax[i][2] ) )


class SOFI(algo):
    """

    The Symmetry Operation FInder (SOFI) is an algorithm for finding point group symmetry operations of
    atomic structures. It solves the problem:

    .. math::
        A = \\theta A

    where :math:`A` is an atomic structure, and :math:`\\theta` is a symmetry operation in the form of
    3x3 orthonormal matrix.

    The main result of SOFI is a list of operations :math:`\\theta`. This list can be post-processed to
    obtain the point group name, and other properties associated to each :math:`\\theta`.

    For further details see the publication:

    M Gunde, N Salles, L Grisanti, L Martin-Samos, A Hemeryck:
    J. Chem. Phys. 2024, 161, 6, 062503
    DOI: https://doi.org/10.1063/5.0215689
    arXiv: https://arxiv.org/abs/2408.06131


    In order to use the module, the path to this file should be added to the environment variable PYTHONPATH:

    .. code-block:: bash

        export PYTHONPATH=/your/path/to/IRA/interface


    Example import:

    >>> import ira_mod
    >>> sofi = ira_mod.SOFI( )

    """

    # nmax is imposed by sofi_tools.f90; it gives the maximal number of symmetries
    # if you need more than than then change in sofi_tools.f90 and here, and recompile sofi
    _nmax=400

    def compute(self, nat, typ_in, coords, sym_thr, origin=None, prescreen_ih=False ):
        """
        .. _sofi.compute:

        this is a wrapper to libira_compute_all() from library_sofi.f90
        Description: This routine computes all the relevant PG symmetry data of the input structure.
        Unless an `origin` point is specified, the input structure is shifted to the geometric center by
        this python interface.

        **== input: ==**

        :param nat: number of atoms in the structure
        :type nat: int

        :param typ_in: atomic types of the atoms
        :type typ_in: string or int

        :param coords: the atomic positions of the structure
        :type coords: np.ndarray((nat,3), dtype=float)

        :param sym_thr: symmetry threshold value
        :type sym_thr: float

        **== optional: ==**

        :param origin: if specified, it will be taken as origin point, otherwise the geometric center is the origin
        :type origin: np.ndarray(3, dtype=float)

        :param prescreen_ih: if True, force sofi to check for Ih group already during computation.
        :type prescreen ih: bool

        **== output: ==**

        instance of ``ira_mod.sym_data()`` object, which cotains the full symmetry data of the input structure.
        See `sym_data`_.

        object `sym_data` contains the full symmetry data of the input structure.
        To access info about symmetry element *i*, get the data as: `sym_data.var[i]`, where `var`
        is one of the variables in the `sym_data` object:

            `n_sym`, `matrix`, `perm`, `op`, `n`, `p`, `axis`, `angle`, `dHausdorff`, `pg`, `prin_ax`

        see help( ira_mod.sym_data ) for full description of these variables.


        """

        ## initialize instance of sym_data
        sym = sym_data()

        ## nmax value
        nmax=self._nmax

        # check if input types are already c_int or not
        typ = self.tf_int( typ_in )

        # if origin is not provided, use geometric center
        gc = np.mean( coords, axis=0 )
        # if origin is provided, use that
        if( origin is not None ):
            # check size, type
            o_size = np.size(origin)
            if( o_size != 3 ):
                msg = "the `origin` must be a 3-vector! Received size:", o_size
                raise ValueError( msg )
            if( not isinstance(origin[0], (float, np.float32, np.float64, int, np.int32, np.int64)) ):
                msg = "unexpected datatype of `origin`! Expected: float or int"
                raise ValueError( msg )
            ## set the input origin as gc
            gc = origin

        # print( "gc used:", gc )
        # shift coords to origin
        coords = coords - gc

        # for i in range(nat):
        #     print( typ[i], *coords[i] )

        # input data
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        thr = c_double( sym_thr )
        check_ih = c_bool( prescreen_ih )


        # allocate output data
        nmat = c_int()
        mat_data = (c_double*9*nmax)()
        perm_data = (c_int*nat*nmax)()
        op_data = (c_char*1*nmax)()
        n_data = (c_int*nmax)()
        p_data = (c_int*nmax)()
        ax_data = (c_double*3*nmax)()
        angle_data = (c_double*nmax)()
        dHausdorff_data = (c_double*nmax)()
        pg = (c_char*11)()
        n_pax = c_int()
        pax_data = (c_double*3*nmax)()
        cerr = c_int()


        # have to set argtypes in here, since nat is not know in init
        self.lib.libira_compute_all.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), c_double, c_bool, \
              POINTER(c_int), POINTER(POINTER(c_double*9*nmax)), POINTER(POINTER(c_int*nat*nmax)), \
              POINTER(POINTER(c_char*1*nmax)), POINTER(POINTER(c_int*nmax)), POINTER(POINTER(c_int*nmax)), \
              POINTER(POINTER(c_double*3*nmax)), POINTER(POINTER(c_double*nmax)), POINTER(POINTER(c_double*nmax)), \
              POINTER(POINTER(c_char*11)), POINTER(c_int), POINTER(POINTER(c_double*3*nmax)), POINTER(c_int) ]
        self.lib.libira_compute_all.restype=None

        # call routine from library.f90
        self.lib.libira_compute_all( n, t, c, thr, check_ih, \
                                  nmat, pointer(mat_data), pointer(perm_data), \
                                  pointer(op_data), pointer(n_data), pointer(p_data), \
                                  pointer(ax_data), pointer(angle_data), pointer(dHausdorff_data), pointer(pg), \
                                  pointer(n_pax), pointer(pax_data), cerr )
        if cerr.value != 0:
            raise ValueError("nonzero error value obtained from libira_compute_all()")

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

        sym.op = np.empty(n_op,dtype="U1")
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

        sym.dHausdorff=np.zeros(n_op, dtype=float)
        for n in range(n_op):
            sym.dHausdorff[n]=dHausdorff_data[n]
        # sym.dHausdorff = dHausdorff_data[:n_op]

        sym.n_prin_ax = n_pax.value
        sym.prin_ax=np.zeros((sym.n_prin_ax, 3), dtype=float)
        for n in range(sym.n_prin_ax):
            sym.prin_ax[n]=pax_data[n][:]

        # return the instance of sym_data
        return sym


    def get_symm_ops(self, nat, typ_in, coords, sym_thr, prescreen_ih=False ):
        """
        Wrapper to the libira_get_symm_ops routine from library_sofi.f90.
        Description: finds the symmetry operation maytices of the structure in input,
        which fit the threshold `sym_thr`.

        **== input: ==**

        :param nat: number of atoms in the structure
        :type nat: int

        :param typ_in: atomic types of the atoms
        :type typ_in: string or int

        :param coords: the atomic positions of the structure
        :type coords: np.ndarray((nat,3), dtype=float)

        :param sym_thr: symmetry threshold value
        :type sym_thr: float

        **== optional ==**

        :param prescreen_ih: if True, force sofi to check for Ih group already during computation.
        :type prescreen ih: bool


        **== output: ==**

        :param n_sym: number of symmetry elements found
        :type n_sym: int

        :param matrix: the list of matrices of symmetry operations
        :type matrix: np.ndarray((n_sym, 3, 3), dtype = float)

        """

        nmax=self._nmax

        # check if input types are already c_int or not
        typ = self.tf_int( typ_in )

        # input data
        n = c_int( nat )
        t = typ.ctypes.data_as( POINTER(c_int) )
        c = coords.ctypes.data_as( POINTER(c_double) )
        thr = c_double( sym_thr )
        check_ih = c_bool( prescreen_ih )
        cerr = c_int()


        # output data
        nmat = c_int()
        mat_data = (c_double*9*nmax)()
        # have to set argtypes in here, since nat is not know in init
        self.lib.libira_get_symm_ops.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), c_double, c_bool, \
              POINTER(c_int), POINTER(POINTER(c_double*9*nmax)), POINTER(c_int) ]
        self.lib.libira_get_symm_ops.restype=None

        self.lib.libira_get_symm_ops( n, t, c, thr, check_ih, nmat, pointer(mat_data), cerr )
        if cerr.value != 0:
            raise ValueError( "nonzero error value obtained form libira_get_symm_ops()")

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
        wrapper to libira_get_pg() from library_sofi.f90
        Description: find the PG of input list of operations. Returns also the
        list of all equivalent principal axes.

        **== input: ==**

        :param nm_in: number of matrices in the `mat_list` array
        :type nm_in: integer

        :param mat_list: list of matrices containging symmetry operations
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )

        **== optional ==**

        :param verb: verbosity of the output, verb=True outputs a report from get_pg() routine.
        :type verb: logical


        **== output: ==**

        :param pg: associated Point Group (PG)
        :type pg: string

        :param n_prin_ax: number of equivalent principal axes
        :type n_prin_ax: integer

        :param prin_ax: list of principal axes of the PG
        :type prin_ax: np.ndarray( (n_prin_ax, 3), dtype = float )

        """
        # input data
        n = c_int( nm_in )
        mats = mat_list.ctypes.data_as( POINTER(c_double) )
        cverb = c_bool( verb )
        cerr = c_int()


        # output data
        pg = (c_char*11)()
        cnprin_ax = c_int()
        pprin_ax = (c_double*3*nm_in)()

        # have to set argtypes in here, since nat is not know in init
        self.lib.libira_get_pg.argtypes = \
            [ c_int, POINTER(c_double), POINTER(POINTER(c_char*11)), \
              POINTER(c_int), POINTER(POINTER(c_double*3*nm_in)), \
              c_bool, POINTER(c_int) ]
        self.lib.libira_get_pg.restype=None
        self.lib.libira_get_pg( n, mats, pointer(pg), pointer(cnprin_ax), pointer(pprin_ax), cverb, cerr)
        if cerr.value != 0 :
            raise ValueError("nonzero error value ontained from libira_get_pg()")

        pg=pg.value.decode()

        n_prin_ax = cnprin_ax.value
        prin_ax=np.zeros((n_prin_ax, 3), dtype=float)
        for n in range(n_prin_ax):
            prin_ax[n]=pprin_ax[n][:]

        return pg, n_prin_ax, prin_ax

    def get_unique_ax_angle( self, nm_in, mat_list ):
        """
        wrapper to libira_get_unique_ax_angle() from library_sofi.f90
        Description: input list of symmetry matrices, output lists
        of op, ax, angle for each symmetry operation in the list, such that
        axes are aligned, and angle has a correspoding +/- sign.

        This is a wrapper to libira_unique_ax_angle() routine from library_sofi.f90

        **== input: ==**

        :param nm_in: number of matrices in the `mat_list` array
        :type nm_in: integer

        :param mat_list: list of matrices containging symmetry operations
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )


        **== output: ==**

        :param op: Schoeflies Op name, possible values:
                         - "E" = identity
                         - "I" = inversion
                         - "C" = rotation
                         - "S" = (roto-)reflection
        :type op: array of length-1 strings, nd.ndarray( nm_in, dtype = "U1" )

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
        nl = n_op+1

        # output
        op_data = (c_char*nl)()
        ax_data = (c_double*3*n_op)()
        angle_data = (c_double*n_op)()
        cerr = c_int()

        self.lib.libira_unique_ax_angle.argtypes = \
            [ c_int, POINTER(c_double), POINTER(POINTER(c_char*nl)) , \
              POINTER(POINTER(c_double*3*n_op)), POINTER(POINTER(c_double*n_op)), POINTER(c_int) ]
        self.lib.libira_unique_ax_angle.restype=None

        self.lib.libira_unique_ax_angle( n_m, mats, pointer(op_data), \
                                      pointer(ax_data), pointer(angle_data), pointer(cerr) )
        if cerr.value != 0:
            raise ValueError( "error in unique_ax_angle")

        op = np.empty(n_op,dtype="U1")
        for n in range(n_op):
            op[n] = op_data[1*n:1*n+1].decode()
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
        Description: analyse the input matrix, output the Schoenflies notation: Op n^p, the
        axis, and angle.
        The notation can be transformed into an angle in radians, or degree:

            angle_radian = p / n * 2pi
            angle_degree = p / n * 360

        e.g. C 8^3 corresponds to 3/8th of full circle = 135 degress = 2.3563 radian, p/n = 0.375

        This is a wrapper to libira_analmat() routine from library_sofi.f90


        **== input: ==**

        :param rmat: input matrix
        :type rmat: np.ndarray((3, 3), dtype = float )


        **== output: ==**

        :param op: Schoeflies Op name, possible values:
                         - "E" = identity
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
        :type angle: np.float32

        """

        # input data
        mat = rmat.ctypes.data_as( POINTER(c_double) )

        # output
        op=(c_char*2)()
        n=c_int()
        p=c_int()
        cax=(c_double*3)()
        angle=c_double()
        cerr = c_int()

        self.lib.libira_analmat.argtypes=\
            [ POINTER(c_double), POINTER(POINTER(c_char*2)), \
              POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_double*3)), \
              POINTER(c_double), POINTER(c_int) ]
        self.lib.libira_analmat.restype=None

        self.lib.libira_analmat( mat, pointer(op), pointer(n), pointer(p), \
                              pointer(cax), pointer(angle), cerr )
        if cerr.value != 0:
            raise ValueError( "nonzero error value obtained from libira_analmat()")

        op=op.value.decode()
        n=n.value
        p=p.value
        angle=np.float32(angle.value)
        ax = np.zeros([3], dtype=float)
        for i in range(3):
            ax[i] = cax[i]

        return op, n, p, ax, angle


    def ext_Bfield( self, n_mat, mat_list, b_field ):

        """
        wrapper to libira_ext_Bfield() from library_sofi.f90
        Description: Impose a constraint on the relevant symmetry operations,
        in the form of an external magnetic field as arbitrary vector.

        Implementation of the functionality proposed in:


            Ansgar Pausch, Melanie Gebele, Wim Klopper;
            Molecular point groups and symmetry in external magnetic fields.
            J. Chem. Phys. 28 November 2021; 155 (20): 201101.
            DOI: https://doi.org/10.1063/5.0069859


        The magnetic field acts as an "external constraint" on the symmetry matrices,
        and can thus be seen as a post-processing of the unconstrained (full) PG of a structure.

        **== Input: ==**

        :param n_mat: number of matrices in the input `mat_list`
        :type n_mat: integer

        :param mat_list: list of 3x3 matrices of symmetry operations in input
        :type mat_list: np.ndarray((n_mat, 3, 3), dtype = float )

        :param b_field: direction of magnetic field
        :type b_field: 3d vector, np.ndarray(3, dtype = float )


        **== output: ==**

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

        self.lib.libira_ext_bfield.argtypes = [ c_int, POINTER(c_double), \
                                             POINTER(c_double), POINTER(c_int), \
                                             POINTER(POINTER(c_double*9*n_mat))]
        self.lib.libira_ext_bfield.restype=None

        self.lib.libira_ext_bfield( n_m, mats, bf, pointer(n_out), pointer(m_out) )

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

        This is a wrapper to libira_get_perm() from library_sofi.f90

        **== Input: ==**

        :param nat: number of atoms in the structure
        :type nat: integer

        :param typ_in: atomic types of the structure
        :type typ_in: np.ndarray( nat, dtype = integer or string )

        :param coords: atomic positions in the structure
        :type coords: np.ndarray((nat, 3), dtype = float )

        :param nm_in: number of matrices in the `mat_list` input
        :type nm_in: integer

        :param mat_list: list of input matrices
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )

        **== Output: ==**

        :param perm: list of permutations corresponding to each matrix operation in input
        :type perm: np.ndarray( (nm_in, nat), dtype = integer )

        :param dHausdorff: maximal distance (Hausdorff distance) corresponding to the application off
                     each matrix in the `mat_list` to the structure. Can be seen as "score" of each symmetry
        :type dHausdorff: np.ndarray( nm_in, dtype = float )

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
        dHausdorff_data = (c_double*nm_in)()


        # have to set argtypes in here, since nat is not know in init
        self.lib.libira_get_perm.argtypes = \
            [ c_int, \
              POINTER(c_int), \
              POINTER(c_double), \
              c_int, \
              POINTER(c_double), \
              POINTER(POINTER(c_int*nat*nm_in)), \
              POINTER(POINTER(c_double*nm_in)) ]
        self.lib.libira_get_perm.restype=None

        self.lib.libira_get_perm( n, t, c, n_m, mats, \
                               pointer(perm_data), pointer(dHausdorff_data) )

        # cast the result into readable things
        perm = np.ndarray((nm_in, nat),dtype=int)
        for n in range(nm_in):
            perm[n]=perm_data[n][:]
        dHausdorff=np.ndarray(nm_in, dtype=float)
        for n in range(nm_in):
            dHausdorff[n] = dHausdorff_data[n]

        return perm, dHausdorff

    def get_combos( self, nat, typ_in, coords, nm_in, mat_list ):

        '''
        Description:
        Obtain all unique combinations of input matrices, which are symmetry operations of given structure.

        This is a wrapper to the routine libira_get_combos() from library_sofi.f90

        **== Input: ==**

        :param nat: number of atoms in the structure
        :type nat: integer

        :param typ_in: atomic types of the structure
        :type typ_in: np.ndarray( nat, dtype = integer or string )

        :param coords: atomic positions in the structure
        :type coords: np.ndarray((nat, 3), dtype = float )

        :param nm_in: number of matrices in the `mat_list` input
        :type nm_in: integer

        :param mat_list: list of input matrices
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )


        **== Output: ==**

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
        cerr = c_int()

        self.lib.libira_get_combos.argtypes=\
            [ c_int, \
              POINTER(c_int), \
              POINTER(c_double), \
              c_int, \
              POINTER(c_double), \
              POINTER(c_int), \
              POINTER(POINTER(c_double*9*nmax)), \
              POINTER(c_int) ]
        self.lib.libira_get_combos.restype=None

        self.lib.libira_get_combos( n, t, c, nm, mats, nb_out, pointer(cmat_out), cerr )

        if cerr.value != 0:
            raise ValueError("nonzero error value obtained from libira_get)combos()")

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

        This is a wrapper to libira_try_mat() routine from library_sofi.f90

        **== Input: ==**

        :param nat: number of atoms in the structure
        :type nat: integer

        :param typ_in: atomic types of the structure
        :type typ_in: np.ndarray( nat, dtype = integer or string )

        :param coords: atomic positions in the structure
        :type coords: np.ndarray((nat, 3), dtype = float )

        :param theta: the input matrix of an orthonormal operation to be tested on structure
        :type theta: np.ndarray((3,3), dtype = float)


        **== Output: ==**

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

        self.lib.libira_try_mat.argtypes = \
            [ c_int, POINTER(c_int), POINTER(c_double), \
              POINTER(c_double), POINTER(c_double), POINTER(POINTER(c_int*nat)) ]
        self.lib.libira_try_mat.restype = None

        self.lib.libira_try_mat( n, t, c, r, pointer(c_dh), pointer(c_perm) )


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

        This is a wrapper to libira_construct_operation() routine from library_sofi.f90

        **== Input: ==**

        :param op: Schoeflies Op name, possible values:
                         - "E" = identity
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


        **== Output: ==**

        :param matrix: the matrix representation of desired operation
        :type matrix: 3x3 matrix, np.ndarray((3, 3), dtype = float )


        .. note::
            The operation "E" is equivalent to "C" about any axis for angle=0.0; and
            similarly the operation "I" is equivalent to "S" with angle=0.5 about any axis.
            The operation "S" with angle=0.0 is a reflection.


        '''
        c_op=op.encode()
        c_ax = axis.ctypes.data_as(POINTER(c_double))
        c_angle = c_double(angle)

        mat_out=(c_double*9)()
        cerr = c_int()

        self.lib.libira_construct_operation.argtypes = \
            [ c_char_p, POINTER(c_double), c_double, POINTER(POINTER(c_double*9)), POINTER(c_int) ]
        self.lib.libira_construct_operation.restype = None

        self.lib.libira_construct_operation( c_op, c_ax, c_angle, pointer(mat_out), cerr )
        if cerr.value != 0:
            raise ValueError( "nonzero error value obtained from libira_construct_operation()")

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

        This is a wrapper to routine libira_mat_combos() from library_sofi.f90

        **== Input: ==**

        :param nm_in: number of matrices in the `mat_list` input
        :type nm_in: integer

        :param mat_list: list of input matrices
        :type mat_list: np.ndarray( (nm_in, 3, 3), dtype = float )


        **== Output: ==**

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

        self.lib.libira_mat_combos.argtypes=\
            [ \
              c_int, \
              POINTER(c_double), \
              POINTER(c_int), \
              POINTER(POINTER(c_double*9*nmax)) \
             ]
        self.lib.libira_mat_combos.restype=None

        self.lib.libira_mat_combos( nm, mats, nb_out, pointer(cmat_out) )

        nm_out=nb_out.value
        mat_out=np.ndarray( (nm_out,3,3), dtype=float)
        for n in range(nm_out):
            m=0
            for i in range(3):
                for j in range(3):
                    mat_out[n][i][j] = cmat_out[n][m]
                    m+=1
        return nm_out, mat_out

    def matrix_distance( self, m1, m2 ):
        '''
        Compute the distance between two matrices using the ``matrix_distance`` function.
        This function is used internally in SOFI to determine if two matrices are equal or not,
        if the value of distance is below ``m_thr``, then the matrices are considered equal.
        By default, ``m_thr = 0.73`` which should be enough to distinguish matrices separated by
        operation equivalent to C 12^1. These values can be found in sofi_tools.f90

        **== input ==**

        :param m1: matrix 1
        :type m1: np.ndarray( [3, 3], dtype = float)

        :param m2: matrix 2
        :type m2: np.ndarray( [3, 3], dtype = float)

        **== output ==**

        :param d: distance value
        :type d: float

        '''

        # input
        mc1 = m1.ctypes.data_as(POINTER(c_double))
        mc2 = m2.ctypes.data_as(POINTER(c_double))

        # output
        cd = c_double()

        self.lib.libira_matrix_distance.restype=None
        self.lib.libira_matrix_distance.argtypes=[POINTER(c_double), POINTER(c_double), POINTER(c_double)]

        self.lib.libira_matrix_distance( mc1, mc2, pointer(cd) )

        d = cd.value
        return d


    def check_collinear( self, nat, coords ):
        '''
        Check if the input structure is collinear or not.

        **== input ==**

        :param nat: number of atoms in the structure
        :type nat: int

        :param coords: the atomic positions of the structure
        :type coords: np.ndarray((nat,3), dtype=float)

        **== output ==**

        :param collinear: true if structure is collinear
        :type collinear: bool

        :param ax: axis of collinearity, if structure is collinear
        :type ax: np.ndarray(3, dtype=float)

        '''
        #input
        n = c_int( nat )
        c = coords.ctypes.data_as( POINTER(c_double) )

        # output
        c_col = c_bool()
        c_ax = (c_double*3)()

        self.lib.libira_check_collinear.restype = None
        self.lib.libira_check_collinear.argtypes = [ c_int, POINTER(c_double), POINTER(c_bool), \
                                                     POINTER(POINTER(c_double*3))]

        self.lib.libira_check_collinear( n, c, pointer(c_col), pointer(c_ax) )
        ax = np.zeros([3], dtype=float)
        for i in range(3):
            ax[i] = c_ax[i]
        return c_col.value, ax
