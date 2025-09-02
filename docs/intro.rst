.. _intro:


Introduction
============

Even though this library is called the **IRA library**,
there are three main algorithms contained inside:

 * CShDA: Constrained Shortest Distance Assignments;
 * IRA: Iterative Rotations and Assignments;
 * SOFI: Symmetry Operation FInder.

Each of those algorithms is designed to solve a different problem related to characterization of generic atomic
structures.

Some other miscallaneous routines are contained, see the section of :ref:`How-to guides <sec_ex_howto>`.


CShDA
-----

Constrained Shortest Distance Assignments (CShDA, [ira2021]_) is the really fundamental algorithm for the rest of this library.
It solves the so-called Linear Assignment Problem (LAP), also referred to as the "point-set registry",
and then computes the distance between two atomic structures :math:`A` and :math:`B`.

More precisely, it solves the problem sometimes referred to as `bottleneck-LAP`, and uses that solution to
compute the Hausdorff distance:

.. math::
   d_H( A, B ) = \max_i \min_j \big\{ d(a_i, b_j) \big\} \qquad(1)

where :math:`a_i \in A` and :math:`b_j \in B` are atomic positions, and :math:`d()` is the Cartesian
distance function.

By first solving the LAP, computing the distance between two atomic structures becomes trivial.
In fact, such distance is invariant to permutations of atoms :math:`P_B`, and variant to all other
rigid transformations.
CShDA returns the atomic permutations, and values of atom-atom distances w.r.t. the permutation:
distances :math:`d(a_i, b_i)`,
with indices :math:`i`, after applying the permutation computed in the beginning (LAP solution).
Finally, the distance between :math:`A` and :math:`B` is taken as the maximal value in the set.

It enforces a strict one-to-one assignment of the atoms, and works also for structures containing
different numbers of atoms.

Since CShDA can be relatively heavy to compute many times, a heuristic early-exit criterion is set up, and used whenever possible.
The criterion is in the form of a distance threshold: as soon as it is established that the final distance
value returned by CShDA cannot be below the given threshold, the computation exits with a high distance value.
The heuristic can be disregarded by simply inputting a very high value for this threshold.

(under construction)


IRA
---

Iterative Rotations and Assignments (IRA, [ira2021]_) solves the so-called `shape matching` problem, also sometimes called `structural alignment`, or similar:

.. _eq-pb3:
.. math::
   P_B B = \mathbf{R} A + \mathbf{t} \qquad (2)

where :math:`A` and :math:`B` are two atomic structures, :math:`\mathbf{R}` is a rotation/reflection matrix,
:math:`\mathbf{t}` is a translation vector, and :math:`P_B` is a permutation matrix of atomic indices in :math:`B`.

The problem of finding an optimal rotation :math:`\mathbf{R}` when :math:`A` and :math:`B` do *not* include permutations is known as
the `orthogonal Procrustes problem`, and can be easily solved by Singular Value Decomposition method (SVD),
see `the wikipedia article <https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem>`_.
The routines to solve it are also included in this library.

The shape-matching problem of Eq. :ref:`(2) <eq-pb3>` is however more complicated since :math:`A` and :math:`B` generally *do* include permutations,
but it can be rewritten as an optimization problem:

.. _eq-argmin:
.. math::
   \DeclareMathOperator*{\argmin}{arg\,min}
   \argmin_{\bf R,t} \big\{ D({\bf R}A + {\bf t}, B) \big\} \qquad (3)

in which :math:`D` is a general distance function between two sets, that is (i) variant under :math:`{\bf R}` and
:math:`{\bf t}`, (ii) invariant under permutation :math:`P_B`, and (iii) returns value 0 when :math:`{\bf R}` and
:math:`{\bf t}` are such that Eq. :ref:`(2) <eq-pb3>` is satisfied, i.e. when the best match is found.

The IRA algorithm solves the Eq. :ref:`(3) <eq-argmin>` by parsing the space of rigid transformations associated to the atomic structure,
computing :math:`D()` for each candidate, and returning the single transformation that minimizes it.
In presence of distortions, it calls the SVD-based algorithm to further minimize the rotations.
The function :math:`D()` in Eq. :ref:`(3) <eq-argmin>` is computed by CShDA.

When :math:`A` and :math:`B` are congruent, IRA is guaranteed to return the optimal solution,
independently of the initial orientation and permutation of the structures, which is not entirely
obvious (see benchmarks in [ira2021]_).
In case of nearly-congruent structures, this becomes debatable, since
the notion of near-congruency is not uniquely defined, and is dependent on the involved structures.

Therefore, IRA is not recommended for applications where the `value` of (dis-)similarity between structures
is the key quantity. It is however highly efficient in tasks like retrieving a certain structure from a
list of structures, where the structure is certain to be present in the list
(in any orientation and permutation).

.. [ira2021] M. Gunde, N. Salles, A. Hemeryck, L. Martin-Samos, `JCIM, 2021`, `DOI: 10.1021/acs.jcim.1c00567 <https://doi.org/10.1021/acs.jcim.1c00567>`_, `arXiv: 2111.00939 <https://export.arxiv.org/abs/2111.00939>`_



SOFI
----

Symmetry Operation Finder (SOFI, [sofi2024]_) is an algorithm for finding point group symmetry operations of
atomic structures.

By definition, the transformation of a structure by a symmetry operation, should give a structure that
is equivalent to the original, up to a permutation of indices:

.. _eq-sofi2:
.. math::
   P_A A = \boldsymbol{\theta} A \qquad (4)

where :math:`A` is an atomic structure,
:math:`P_A` is a permutation of atomic indices in :math:`A`, and
:math:`\boldsymbol{\theta}` is a symmetry operation in the form of
3x3 orthonormal matrix.

Notice the similarity of Eq. :ref:`(4) <eq-sofi2>` to Eq. :ref:`(2) <eq-pb3>`: the structure :math:`B` is now
equal to :math:`A`, the rigid transformation :math:`\mathbf{R}` becomes a symmetry
operation :math:`\boldsymbol{\theta}`, and the translation :math:`\mathbf{t}` vanishes.

When the structure :math:`A` contains point group symmetries, Eq. :ref:`(4) <eq-sofi2>` has degenrate
solutions in form of pairs :math:`(\boldsymbol{\theta}, P_A)`.
The set of all such pairs represents the set of point group symmetry operations for the structure.
SOFI solves this problem.
It can be seen as an extension of IRA, where IRA gives a single, optimal solution to matching two (near-)congruent
structures, SOFI gives all degenerate solutions of matching a structure to itself.

.. [sofi2024] M. Gunde, N. Salles, L. Grisanti, L. Martin-Samos, A. Hemeryck, `JCP, 2024`, `DOI: 10.1063/5.0215689 <https://doi.org/10.1063/5.0215689>`_, `arXiv: 2408.06131 <https://arxiv.org/abs/2408.06131>`_

