.. _sofi_howto:

###############################
SOFI tutorial and How-to guides
###############################


.. note::

   The guides are given for python use. The same can be done from C/Fortran by calls to
   appropriate routines, which generally correspond in the name to the ones from python.
   For more info refer to :ref:`ref_api` and :ref:`src_refs`.




.. contents:: Contents
   :local:
   :depth: 2


Importing the ira module
========================

The IRA library is imported into python by:

>>> import ira_mod


.. warning::
   If the ``ira_mod`` module cannot be found at ``import``, then make sure there is a path to ``/IRA_library/interface``
   in the environment variable ``PYTHONPATH``.

   .. code-block:: bash

      echo $PYTHONPATH

   If not, add it by:

   .. code-block:: bash

      export PYTHONPATH=/your/path/to/IRA_library/interface:$PYTHONPATH

   Alternatively, use ``pip`` to install the module:

   .. code-block:: bash

      python -m pip install .


The corresponding algorithm class (IRA or SOFI) has to be initialised by:

>>> ira = ira_mod.IRA()

or

>>> sofi = ira_mod.SOFI()

now, the functions in either class are availably by typing ``ira.<funtion_name>`` or ``sofi.<function_name>``.
Quick help can be accessed by ``help( ira )``, ``help( ira.<function_name> )`` or the same for ``sofi``.


Construct matrices from axis-angle representation
=================================================

In the Schoenflies notation, orthonormal transformation matrices can be written in the generic form ``Op n^p``,
where the ``Op`` specifies the `action` of a matrix,
and ``n^p`` specify the angle :math:`\alpha` as: :math:`{p}/{n} = \alpha`.
Possible values for ``Op`` are:

 - ``E`` for identity,
 - ``C`` for proper rotation,
 - ``S`` for reflection (with zero angle), or improper rotation (with non-zero angle),
 - ``I`` for point-inversion.

To generate a matrix from its Schoenflies notation (and an axis), for example S 12^5 along axis (0, 0, 1), the function ``construct_operation()`` can be used, specifying the axis ``ax`` by:

   >>> import numpy as np
   >>> ##
   >>> ## define the axis
   >>> ax = np.array([ 0., 0., 1.] )
   >>> ##
   >>> ## construct operation S 12^5 along that axis, with angle 5/12
   >>> matrix = sofi.construct_operation( "S", ax, 5/12 )
   >>> matrix
   array([[-0.8660254, -0.5      ,  0.       ],
          [ 0.5      , -0.8660254,  0.       ],
          [ 0.       ,  0.       , -1.       ]])
   >>> ##
   >>> ## giving negative angle should return the transposed matrix
   >>> sofi.construct_operation( "S", ax, -5/12 )
   array([[-0.8660254,  0.5      ,  0.       ],
          [-0.5      , -0.8660254,  0.       ],
          [ 0.       ,  0.       , -1.       ]])
   >>> ##
   >>> ## giving negative axis should return the transpose also
   >>> sofi.construct_operation( "S", -ax, 5/12 )
   array([[-0.8660254,  0.5      ,  0.       ],
          [-0.5      , -0.8660254,  0.       ],
          [ 0.       ,  0.       , -1.       ]])


.. note::
   The axis ``ax`` on input does not need to be normalised.


.. _analmat:

Obtain axis-angle and Schoenflies fom matrix
============================================

An orthonormal 3x3 matrix can be analysed to obtain its Schoeflies representation of the format ``Op n^p``,
and the axis-angle representation by calling the ``analmat()`` function:

   >>> ## create a matrix for C 5^2 along axis (1., -1., 1.)
   >>> matrix = sofi.construct_operation( "C", np.array([1., -1., 1.]), 2/5 )
   >>> ##
   >>> ## analyse it
   >>> sofi.analmat( matrix )
   ('C', 5, 2, array([ 0.57735027, -0.57735027,  0.57735027]), 0.4)
   >>> ## save the output
   >>> op, n, p, ax, angle = sofi.analmat( matrix )

The Schoeflies symbol is then ``Op n^p``. The ``angle`` is in units of :math:`2\pi`, i.e. ``angle=0.5`` is half
the full circle. The axis ``ax`` on output is normalised.

.. note::
   the axis ``ax`` comes from a diagonalisation procedure, therefore any :math:`\pm` direction is a
   valid solution. To remove this ambiguity, the convention is that the axis is flipped such that its components are
   :math:`z>0`, if :math:`z=0` then :math:`x>0`, and if :math:`x=0` then :math:`y>0` (all within
   threshold of numerical precision, which is ``epsilon=1e-6`` by default). The orientation of the angle is then decided based on this axis convention.
   Therefore it can happen that analysis of a matrix constructed as:
   
      >>> matrix = sofi.construct_operation( "C", np.array([-0.3, 1., 0.]), 3/8 )
   
   will flip its axis and angle :
   
      >>> sofi.analmat( matrix )
      ('C', 8, 3, array([ 2.87347886e-01, -9.57826285e-01, -1.60749682e-16]), -0.375 )

.. warning::
   The computation of ``n`` and ``p`` in SOFI is limited to a certain order, which is by default 200 at maximum.
   If the order of a matrix is larger than that, ``analmat`` will return ``n`` and ``p`` which are wrong, but
   as close as possible to truth, within the `resolution` of 1/200. The ``angle`` will have
   the correct value in any case.
   In order to modify this behaviour, edit the ``lim_n_val`` parameter, as described :ref:`here <modif_m_thr>`.


Generate a cyclic group from combinations of matrices
=====================================================

Two or more matrices can be used to create a cyclic group. A cyclic group means any combination of the
elements always generates an element that is inside the group. This can be done by calling the ``mat_combos()``
function:

.. code-block:: python

   >>> ## create an empty list of two 3x3 matrices
   >>> mat_list = np.zeros( [2, 3, 3], dtype=float)
   >>> ##
   >>> ## the first matrix flips over x, and the second over z
   >>> mat_list = np.array([[[-1.,  0.,  0.],
   ...                       [ 0.,  1.,  0.],
   ...                       [ 0.,  0.,  1.]],
   ...
   ...                       [[ 1.,  0.,  0.],
   ...                        [ 0.,  1.,  0.],
   ...                        [ 0.,  0., -1.]]])
   >>> ##
   >>> ## create combinations until group completeness
   >>> n_combo, combo_list = sofi.mat_combos( 2, mat_list )
   >>> n_combo
   4
   >>> combo_list
   array([[[-1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]],
           
          [[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0., -1.]],

          [[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]],

          [[-1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0., -1.]]])



Determine point group from a list of matrices
=============================================

A point group can be deduced from list of 3x3 orthonormal matrices, using the ``get_pg()`` function.
The determination follows the standard flowchart, i.e. https://symotter.org/assets/flowchart.pdf

   >>> ## create an empty list of four 3x3 matrices
   >>> mat_list = np.zeros( [4, 3, 3], dtype=float)
   >>> ##
   >>> ## add some operations:
   >>> ## identity
   >>> mat_list[0] = sofi.construct_operation("E", np.array([1., 0., 0.]), 0)
   >>> ## mirror over x
   >>> mat_list[1] = sofi.construct_operation("S", np.array([1., 0., 0.]), 0)
   >>> ## mirror over y
   >>> mat_list[2] = sofi.construct_operation("S", np.array([0., 1., 0.]), 0)
   >>> ## mirror over z
   >>> mat_list[3] = sofi.construct_operation("S", np.array([0., 0., 1.]), 0)
   >>> ##
   >>> ## create complete cyclic group by combinations
   >>> n_combo, combo_list = sofi.mat_combos( 4, mat_list )
   >>> ##
   >>> ## what operations does the new list contain?
   >>> for mat in combo_list:
   ...    sofi.analmat( mat )
   ... 
   ('E', 0, 1, array([1., 0., 0.]), 0.0)
   ('S', 0, 1, array([1., 0., 0.]), 0.0)
   ('S', 0, 1, array([0., 1., 0.]), 0.0)
   ('S', 0, 1, array([0., 0., 1.]), 0.0)
   ('C', 2, 1, array([0., 0., 1.]), 0.5)
   ('C', 2, 1, array([0., 1., 0.]), 0.5)
   ('C', 2, 1, array([1., 0., 0.]), 0.5)
   ('I', 2, 1, array([1., 0., 0.]), 0.5)
   >>> ##
   >>> ## get point group and list of equivalent principal axes of the new list
   >>> pg, n_prin_ax, prin_ax = sofi.get_pg( n_combo, combo_list )
   >>> pg
   'D2h'
   >>> prin_ax
   array([[0., 0., 1.],
          [0., 1., 0.],
          [1., 0., 0.]])
   >>> ##
   >>> ## a more verbose output can be obtained by setting `verb=True`:
   >>> sofi.get_pg( n_combo, combo_list, verb = True )


Test generator elements
=======================

Now we can test by trial-and-error if certain symmetry elements are generator elements of a group.
For example, the Td point group should be possible to generate from two S4 operations on perpendicular axes.

   >>> ## create empty list of two 3x3 matrices
   >>> mat_list = np.zeros( [2, 3, 3] )
   >>> ##
   >>> ## create two S4 operations, on perpendicular axes
   >>> mat_list[0] = sofi.construct_operation("S", np.array([1., 0., 0.]), 1/4)
   >>> mat_list[1] = sofi.construct_operation("S", np.array([0., 1., 0.]), 1/4)
   >>> ##
   >>> ## generate all combinations
   >>> nc, mc = sofi.mat_combos(2, mat_list)
   >>> ##
   >>> ## determine point group
   >>> sofi.get_pg( nc, mc )
   ('Td', 4, array([[-0.57735027, -0.57735027,  0.57735027],
          [ 0.57735027,  0.57735027,  0.57735027],
          [-0.57735027,  0.57735027,  0.57735027],
          [ 0.57735027, -0.57735027,  0.57735027]]))


.. _mat_dist:

Matrix distance, or the resolving power of SOFI
===============================================

In SOFI, two matrices are considered equal when the function ``matrix_distance()`` returns a
value below the threshold ``m_thr``, the default value for which is ``m_thr=0.044``. Example:

   >>> ## create two matrices: S4 and C2 on the same axis
   >>> m1 = sofi.construct_operation( "S", np.array([ 1., 0., 0.]), 1/4 )
   >>> m2 = sofi.construct_operation( "C", np.array([ 1., 0., 0.]), 1/2 )
   >>> ##
   >>> ## compute distance between them
   >>> sofi.matrix_distance( m1, m2 )
   2.8284271247461903

The value of ``matrix_distance`` can be seen as the order of a matrix ``R`` needed to transform ``m1`` into ``m2``.

   >>> ## generate matrix R which transforms m1 into m2:
   >>> R = np.matmul( m1.T, m2 )
   >>> ##
   >>> ## analyse R
   >>> sofi.analmat( R )
   ('S', 4, 1, array([1., 0., 0.]), 0.25)

The threshold ``m_thr`` specifies the maximal order of transformation matrix ``R``, through the computation of the ``matrix_distance()``.
When the distance between two matrices ``m1`` and ``m2`` is above the ``m_thr`` threshold, SOFI will consider the two matrices as different, and when the distance is below ``m_thr``, the matrices are regarded as equal.

This can be seen by constructing two very similar matrices ``m1`` and ``m2``, and computing the matrix ``R`` which transforms one into the other. Thus, ``R`` should be very similar to the identity matrix.
If the analysis of ``R`` returns the identity matrix, then matrices ``m1`` and ``m2`` are considered equal.

   >>> ## create matrices which are similar:
   >>> m1 = sofi.construct_operation( "C", np.array([1., 0., 0.]), 0.5 )
   >>> m2 = sofi.construct_operation( "C", np.array([1., 0., 0.]), 0.503 )
   >>> ## get R
   >>> R = np.matmul( m1.T, m2 )
   >>> sofi.analmat( R )
   ('C', 1, 1, array([1., 0., 0.]), 0.003)
   >>> ## notice C 1^1 is an identity matrix, even if the angle value is in principle correct
   >>> #
   >>> ## compute the distance from m1 to m2
   >>> sofi.matrix_distance( m1, m2 )
   0.026656902985230164
   >>> ## the value is below m_thr=0.044, matrices m1 and m2 are seen as equal



.. note::
   The value of ``m_thr`` effectively determines the `maximal resolving power` of SOFI.
   In case a structure contains symmetry operations with order higher than C200, SOFI will not be able to distinguish them by default.
   If you suspect that is the case, the value of ``m_thr`` can be adjusted to accommodate higher orders, however the ``src`` needs to be recompiled.
   In that case, take care of array sizes, as they might exceed ``nmax``, and to adjust ``lim_n_val``.
   Refer :ref:`here <modif_m_thr>` for more info.




Symmetry operations of an atomic structure
==========================================

Using the ``get_symm_ops()`` function of SOFI to obtain the list of symmetry operations
of a given atomic structure works like:

   >>> import numpy as np
   >>> import ira_mod
   >>> sofi=ira_mod.SOFI()
   >>> ##
   >>> ## create a hypothetical atomic structure with 6 atoms:
   >>> nat = 6
   >>> ## all atomic types equal, integer value 1
   >>> typ = np.ones( [nat], dtype=int)
   >>> ## atomic positions
   >>> coords = np.array([[-0.65 ,  1.126,  0.   ],
   ...                    [-0.65 , -1.126,  0.   ],
   ...                    [ 1.3  , -0.   ,  0.   ],
   ...                    [-1.04 ,  0.   ,  0.   ],
   ...                    [ 0.52 , -0.901,  0.   ],
   ...                    [ 0.52 ,  0.901,  0.   ]])
   >>> ##
   >>> ## specify the symmetry threshold value
   >>> sym_thr = 0.05
   >>> ##
   >>> ## get the symmetry operations in form of 3x3 matrices
   >>> n_mat, mat_list = sofi.get_symm_ops( nat, typ, coords, sym_thr )

The list of matrices can now be input into ``get_pg()``:

   >>> sofi.get_pg( n_mat, mat_list )
   ('D3h', 1, array([[0., 0., 1.]]))

Thus, the structure has D3h point group, with principal axis in the (0, 0, 1) direction.
You can view the hypothetical structure in your favourite visualiser software, and confirm the
symmetry operations and their axes, listed by SOFI:

   >>> for mat in mat_list:
   ...   sofi.analmat( mat )

.. note::
   The structure we have set up as ``coords`` has a geometric mean at (0, 0, 0), it can be confirmed:

      >>> np.mean( coords, axis=0 )
      array([0., 0., 0.])

   In subsequent how-to's we will work with structures where this is not necessarily the case.



Applying symmetry operations
============================

Upon transforming a structure with its symmetry operation, we obtain back the same structure.
Take the same hypothetical structure from before, it has a C3 operation on axis (0, 0, 1):

   >>> ## create a hypothetical atomic structure with 6 atoms:
   >>> nat = 6
   >>> ## all atomic types equal, integer value 1
   >>> typ = np.ones( [nat], dtype=int)
   >>> ## atomic positions
   >>> coords = np.array([[-0.65 ,  1.126,  0.   ],
   ...                    [-0.65 , -1.126,  0.   ],
   ...                    [ 1.3  , -0.   ,  0.   ],
   ...                    [-1.04 ,  0.   ,  0.   ],
   ...                    [ 0.52 , -0.901,  0.   ],
   ...                    [ 0.52 ,  0.901,  0.   ]])
   >>> ##
   >>> ## create C3 along (0, 0, 1)
   >>> c3mat = sofi.construct_operation( "C", np.array([0., 0., 1.]), 1/3)
   >>> ##
   >>> ## create the transformed coords
   >>> coords_tf = np.zeros([nat, 3], dtype=float)
   >>> ##
   >>> ## apply C3 to original coords through np.matmul()
   >>> for i, v in enumerate( coords ):
   ...    coords_tf[i] = np.matmul( c3mat, v )
   ...
   >>> ##
   >>> ## print the transformed structure:
   >>> coords_tf
   array([[-6.504e-01, -1.126e+00,  0.000e+00],
          [ 1.300e+00,  3.576e-05,  0.000e+00],
          [-6.499e-01,  1.126e+00,  0.000e+00],
          [ 5.200e-01, -9.009e-01,  0.000e+00],
          [ 5.205e-01,  9.009e-01,  0.000e+00],
          [-1.040e+00,  7.153e-06,  0.000e+00]], dtype=float)
   >>> ##
   >>> ## notice the vectors are equal (within precision) to the original coords, except permuted.

To obtain the permutation of atoms which happens upon the transformation by a symmetry operation,
SOFI has the ``try_mat()`` function, which returns the value of Hausdorff distance between the original structure,
and the structure transformed by a given matrix, and the corresponding permutation of indices:

   >>> dHausdorff, perm = sofi.try_mat( nat, typ, coords, c3mat )
   >>> ##
   >>> ## print the permutation
   >>> perm
   array([2, 0, 1, 5, 3, 4])
   >>> ## print the Hausdorff distance
   >>> dHausdorff
   0.00033364459005079844


The low value of ``dHausdorff`` confirms that ``c3mat`` is indeed a symmetry operation of the structure defined above.
If you now take ``coords_tf`` from above, permute them by ``perm``, and compute the maximal distance between atoms
``coords[i]`` and ``coords_tf_perm[i]``, you should obtain the value ``dHausdorff``.

   >>> ## permute coords_tf by perm
   >>> coords_tf_perm = coords_tf[ perm ]
   >>> ##
   >>> ## create array for atom-atom distances
   >>> d=np.zeros([nat], dtype=float)
   >>> ##
   >>> ## compute atom-atom distances between the original coords and coords_tf_perm
   >>> for i, v in enumerate( coords ):
   ...    d[i] = np.linalg.norm( v - coords_tf_perm[i] )
   ...
   >>> np.max( d )
   0.000333580064184048


.. note::
   The ``sym_thr`` argument when computing ``get_symm_ops()`` is a threshold in terms of the distance
   ``dHausdorff`` as computed in this section. If an operation returns a distance value beyond ``sym_thr``,
   then SOFI will not consider that operation as a symmetry operation.



Disordered structures: missing operations
=========================================

In case of atomic structures with distortions present in the positions, there could be
some symmetry elements which are either `broken`, or return a distortion higher than expected.
In these cases, SOFI can detect that the number of found symmetry operations does not match
the expected number of operations of the designated point group. The situation can then be resolved
by performing combinations of the found operations, until group completeness.

Set up an atomic structure with distorted atomic positions:

    >>> nat = 21
    >>> typ = np.array([2, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1], dtype=int)
    >>> coords =  np.array([[-0.09854286,  0.07144762, -0.9695    ],
    ...                     [-0.03734286, -1.95445238,  0.7135    ],
    ...                     [-0.00504286, -1.88935238, -1.2304    ],
    ...                     [ 0.02215714, -0.06685238,  1.228     ],
    ...                     [ 1.64625714,  0.96894762, -1.1187    ],
    ...                     [ 1.70545714,  0.90644762,  0.8344    ],
    ...                     [-1.83834286,  1.06694762, -1.1234    ],
    ...                     [-1.67844286,  0.92564762,  0.8333    ],
    ...                     [ 1.74115714, -2.34815238,  1.3447    ],
    ...                     [-1.61704286, -2.87785238,  1.3832    ],
    ...                     [ 1.61885714, -2.79595238, -1.8355    ],
    ...                     [-1.61804286, -2.75785238, -1.8243    ],
    ...                     [ 0.02115714, -0.05535238,  3.2638    ],
    ...                     [ 1.67555714,  2.78904762, -1.7856    ],
    ...                     [ 3.13355714, -0.08455238, -1.7534    ],
    ...                     [ 3.30885714, -0.01745238,  1.4093    ],
    ...                     [ 1.53865714,  2.70804762,  1.4813    ],
    ...                     [-1.51324286,  2.80494762, -1.9041    ],
    ...                     [-3.19054286, -0.07205238, -1.8623    ],
    ...                     [-1.60244286,  2.74904762,  1.4594    ],
    ...                     [-3.21264286, -0.07065238,  1.4563    ]], dtype=float)

View the structure in your visualizer, it should be easy to notice straight away that the (0, 0, 1) axis
should be a C3 axis, however the atomic distortions are relatively large.
Let's set a relatively high symmetry threshold, and try to find the symmetry operations:

   >>> sym_thr = 0.5
   >>> n_mat, mat_list = sofi.get_symm_ops( nat, typ, coords, sym_thr )
   >>> n_mat
   4
   >>> sofi.get_pg( n_mat, mat_list )
   ('C3v-', 1, array([[ 8.61320772e-04, -9.09124124e-03,  9.99958303e-01]]))
   >>> ##

Notice the PG output is ``c3v-``, the minus is a signal that the group
could be identified from the flowchart, but the number of associated
symmetry operations is different than expected for that group. More precisely, the minus sign
indicates that the number is lower than expected. On the contrary, a plus sign would indicate
that SOFI deduced some group, but the number of symmetry elements is higher than expected.

We can now use the ``get_combos()`` function on the list of found symmetries, to form
a complete group of elements that are symmetry elements of atomic structure:

   >>> n_combo, mat_combo = sofi.get_combos( nat, typ, coords, n_mat, mat_list )
   >>> n_combo
   6
   >>> ## two new elements have been generated by combinations. Compute the new PG.
   >>> sofi.get_pg( n_combo, mat_combo )
   ('C3v', 1, array([ 8.61320772e-04, -9.09124124e-03,  9.99958303e-01]))
   >>> ##
   >>> ## the full group has been generated, let's compute permutations and distances
   >>> perm, dHausdorff = sofi.get_perm( nat, typ, coords, n_combo, mat_combo )
   >>> dHausdorff
   array([1.49097439e-15, 4.38744637e-01, 4.35565047e-01, 4.35565047e-01,
          5.12566013e-01, 5.20405469e-01])
   >>> ##
   >>> ## notice the first 4 values are below 0.5 (the sym_thr value used in get_symm_ops),
   >>> ## and the last two which were generated by combinations have ``dHausdorff > 0.5``


And thus we have generated the missing symmetry operations, by performing combinations of the known elements
until group completeness.
The missing operations were not found by SOFI, since their ``dHausdorff`` values are beyond the
``sym_thr=0.5`` we have used in ``get_symm_ops()``, and thus SOFI disregarded them as symmetry elements.

If we repeat the above calculation with ``sym_thr=0.6``, the whole ``C3v`` group should be found straight away.

   >>> sym_thr = 0.6
   >>> n_mat, mat_list = sofi.get_symm_ops( nat, typ, coords, sym_thr )
   >>> sofi.get_pg( n_mat, mat_list )
   ('C3v', 1, array([ 8.61320772e-04, -9.09124124e-03,  9.99958303e-01]))

The feature of performing combinations of elements of a list of matrices gives some flexibility when dealing with
structures with disordered positions, and we do not know the precise value for ``sym_thr`` in advance.


Put it all together: `sofi.compute`
===================================

In order to perform all SOFI computations in one function, that is:
``get_symm_ops()``, then ``get_mat_combos()``, ``get_perm()``, ``analmat()`` and finally ``get_pg()``,
we can simply call the ``compute()`` function:

   >>> sym = sofi.compute( nat, typ, coords, sym_thr )
   >>> ##
   >>> ## see what is in `sym` (use tab)
   >>> sym.
   sym.angle       sym.matrix      sym.n_sym       sym.perm        sym.print()     
   sym.axis        sym.n           sym.op          sym.pg          
   sym.dHausdorff  sym.n_prin_ax   sym.p           sym.prin_ax


The ``compute()`` function returns a ``sym`` object that contains all data computed by SOFI.


Choosing the origin point
=========================

SOFI is agnostic to the choice of the origin point. That means the choice is left to the
user, or application, which calls SOFI.

The most general choice should be the geometric center (arithmetic mean) of the structure, since it is
guaranteed to remain a fixed point for all symmetry elements of the PG of the structure.
The function ``sofi.compute()`` takes an optional argument ``origin``.
If ``origin`` is not specified, SOFI will shift the structure to its geometric center by default.

   >>> ## origin not specified, geo. center will be computed and the structure shifted
   >>> sym = sofi.compute( nat, typ, coords, sym_thr )

In some cases, there can be points other than geometric center, which remain fixed for a subset of the symmetry
elements. These points are then the origin points for subgroups associated to the structure.

If ``origin`` is specified as 3-vector, then that point will be taken as the origin.

Imagine an application where symmetry operations about a given atom are sought, instead of all possible symmetries.
This can be achieved by specifying that point as the ``origin`` in the call to ``sofi.compute()``:

   >>> idx_atm = 7
   >>> my_origin = coords[ idx_atm ]
   >>> ## specify the origin point
   >>> sym = sofi.compute( nat, typ, coords, sym_thr, origin = my_origin )

.. warning::

   The optional argument ``origin`` in ``sofi.compute()`` is available **only in the python interface**,
   the calls to ``sofi_compute_all()`` from other languages do not contain it, meaning you need to shift the structure manually before the call!


HOW-TO: Dealing with linear structures
======================================

Linear structures can have either :math:`C_{\infty v}` or :math:`D_{\infty h}` point groups. The main difference between them is that :math:`D_{\infty h}` has the inversion as symmetry operation, while :math:`C_{\infty v}` does not. The axis of the structure is a rotational axis of infinite order for both groups.

Due to the way the main algorithm of SOFI works, it is limited to structures containing at least 3 noncollinear atoms. Thus, linear structures cannot be explicitly treated with it. The only symmetry operations returned by SOFI when inputting a linear structure will be the identity matrix, and when applicable, the inversion, and reflection over the plane of the axis.


For example, if we create a linear structure without the mirror symmetry, thus group :math:`C_{\infty v}`, SOFI will only find the identity matrix, and the group will be "C1":

.. code-block:: python

   >>> ## create a linear structure with 3 atoms on the x-axis, centered at zero
   >>> nat = 3
   >>> coords = np.array([[-1.0, 0.0, 0.0],
   ...                    [0.0, 0.0, 0.0],
   ...                    [1.0, 0.0, 0.0]])
   >>> ##
   >>> ## specify one of the side atoms as different atomic type
   >>> typ = np.array([1, 1, 2], dtype=int)
   >>> ##
   >>> ## call sofi.compute
   >>> sym = sofi.compute( nat, typ, coords, 0.1 )
   >>> ##
   >>> ## list of matrices has only identity
   >>> sym.matrix
   array([[[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]]])


On the other hand, if we create a structure with inversion and reflection, group :math:`D_{\infty h}`, and call ``sofi.compute()``, the list of matrices has more than one element.

.. code-block:: python

   >>> ## create a linear structure with 4 atoms on the x-axis, already centered at zero
   >>> nat = 4
   >>> coords = np.array([[-1.5,  0. ,  0. ],
   ...                    [-0.5,  0. ,  0. ],
   ...                    [ 0.5,  0. ,  0. ],
   ...                    [ 1.5,  0. ,  0. ]])
   >>> ##
   >>> ## all atoms of the same type:
   >>> typ = np.array([ 1, 1, 1, 1], dtype=int)
   >>> ##
   >>> ## call sofi.compute
   >>> sym = sofi.compute( nat, typ, coords, 0.1 )
   >>> ##
   >>> ## list of matrices has more than one element; notably the mirror over x-axis (the axis of our structure)
   >>> sym.matrix
   array([[[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]],
   
          [[-1.,  0.,  0.],
           [ 0., -1.,  0.],
           [ 0.,  0., -1.]],
   
          [[-1., -0., -0.],
           [ 0.,  1., -0.],
           [-0.,  0.,  1.]],
   
          [[ 1.,  0.,  0.],
           [-0., -1.,  0.],
           [ 0.,  0., -1.]]])


In order to distinguish the linear structures from the others, the library contains a function ``check_collinear()``, which can be used as follows:

.. code-block:: python

   >>> is_collinear, axis = sofi.check_collinear( nat, coords )

Thus if the returned variable ``is_collinear=True``, then the structure in ``coords`` is collinear, and vice versa.
The variable ``axis`` contains the axis of the structure, when it is collinear.

This function can be combined with the ``compute()`` function to properly label point groups of linear structures:

.. code-block:: python

   ## call compute()
   sym = sofi.compute( nat, typ, coords, sym_thr )
   ##
   ## check if structure is collinear
   is_collinear, axis = sofi.check_collinear( nat, coords )
   if( is_collinear ):

      ## if number of found symmetry operations == 1: group should be C_inf_v
      if( sym.n_sym == 1 ):
         ## overwrite the point group as desired
         sym.pg = "Cnv"

      ## more than one foun symmetry operation; group should be D_inf_h
      else:
         sym.pg = "Dnh"



.. _modif_m_thr:

HOW-TO: Modifying the resolving power of SOFI
=============================================

The `maximal resolving power` of SOFI is limited.
In order to modify it, the three parameters in the SOFI source: ``m_thr``, ``nmax``, and ``lim_n_val`` should preferrably be modified, and the source re-compiled. The parameters are located in ``sofi_tools.f90``:

.. code-block:: fortran

  ! real, parameter :: m_thr = 0.73     !! C12
  ! real, parameter :: m_thr = 0.49     !! C18
  ! real, parameter :: m_thr = 0.36     !! C24
  real, parameter :: m_thr = 0.044    !! C200
  ! real, parameter :: m_thr = 0.022    !! C400


The ``m_thr`` is a threshold on matrix distances, its value gives the highest order of an operation that will
still be considered as distinct operation in the SOFI main loop. The value ``m_thr = 0.044`` corresponds
to operation C200, as indicated by the comment in the source.
There are some other values proposed, which can be used by simply
uncommenting them. Value for ``m_thr`` corresponding to other orders of symmetry operations can be computed with
the ``matrix_distance`` function. See also :ref:`here <mat_dist>`.

.. code-block:: fortran

  integer, parameter :: nmax = 400

The ``nmax`` specifies the expected size of input arrays for SOFI. If the number of found symmetry operations
is beyond ``nmax``, SOFI will return an error. Thus if you expect your structure will contain more
than ``nmax`` symmetries, you should edit this value.
Keep in mind that the actual sizes of arrays in the caller software need to be consistent with ``nmax``.

.. code-block:: fortran

  integer, parameter :: lim_n_val = 200

The ``lim_n_val`` is used to find values of ``n`` and ``p`` in ``sofi_analmat()``. If a symmetry operation
with higher order is input, values of ``n`` and ``p`` will be wrong. See also :ref:`here <analmat>`.
