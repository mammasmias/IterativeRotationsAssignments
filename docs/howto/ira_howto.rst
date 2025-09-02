.. _ira_howto:

######################################
IRA & CShDA tutorial and How-to guides
######################################


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


Direct comparison of structures using CShDA
===========================================

Two structures that are already in the same rotational frame, can be directly compared by
computing the atom-atom distances in a way that is invariant to the
permutations with the CShDA algorithm:

   >>> import numpy as np
   >>> ##
   >>> ## set up the first atomic structure containing 10 atoms
   >>> nat1 = 10
   >>> typ1 = np.ones( nat1, dtype=int )
   >>> coords1 = np.array([[-0.49580341,  0.9708181 ,  0.37341428],
   ...                     [-1.05611656, -0.4724503 , -0.37449784],
   ...                     [-0.63509644, -0.66670776,  0.66219897],
   ...                     [-0.83642178,  0.59155936, -0.64507703],
   ...                     [ 0.59636159,  0.80558701,  0.23843962],
   ...                     [ 0.25975284,  0.71540297, -0.78971024],
   ...                     [-0.09743308, -1.03812804, -0.31233049],
   ...                     [ 0.09254502,  0.20016738,  1.03021068],
   ...                     [-0.18424967, -0.24756757, -1.07217522],
   ...                     [ 0.46705991, -0.73516435,  0.56288325]])
   >>> ##
   >>> ## set the second atomic structure as identical, with slightly perturbed atomic positions
   >>> nat2 = nat1
   >>> typ2 = typ1
   >>> coords2= np.array([[-0.50010644,  0.96625779,  0.37944221],
   ...                    [-1.05658467, -0.46953529, -0.37456054],
   ...                    [-0.63373056, -0.66591152,  0.66168751],
   ...                    [-0.83286912,  0.5942803 , -0.64527646],
   ...                    [ 0.59310547,  0.80745772,  0.23711422],
   ...                    [ 0.2636203 ,  0.7126221 , -0.79370807],
   ...                    [-0.09940056, -1.03859144, -0.31064337],
   ...                    [ 0.09208454,  0.19985156,  1.03003579],
   ...                    [-0.18468815, -0.24935304, -1.07257697],
   ...                    [ 0.4691676 , -0.73356138,  0.56184166]])
   >>> ##
   >>> ## add some random permutation to the atoms in second structure
   >>> coords2 = coords2[ [2,4,3,5,9,8,6,7,0,1] ]
   >>> ##
   >>> ## call cshda:
   >>> perm, dist = ira.cshda( nat1, typ1, coords1, nat2, typ2, coords2 )
   >>> ##
   >>> ## the `perm` contains permutations of the second structure, which matches the first structure.
   >>> ## Therefore, the following command should return a structure exactly equal to the first structure:
   >>> coords2[ perm ]
   array([[-0.50010644,  0.96625779,  0.37944221],
          [-1.05658467, -0.46953529, -0.37456054],
          [-0.63373056, -0.66591152,  0.66168751],
          [-0.83286912,  0.5942803 , -0.64527646],
          [ 0.59310547,  0.80745772,  0.23711422],
          [ 0.2636203 ,  0.7126221 , -0.79370807],
          [-0.09940056, -1.03859144, -0.31064337],
          [ 0.09208454,  0.19985156,  1.03003579],
          [-0.18468815, -0.24935304, -1.07257697],
          [ 0.4691676 , -0.73356138,  0.56184166]])
   >>> ##
   >>> ## The `dist` array contains atom-atom distances, upon permuting coords2 by `perm`, such that:
   >>> ## dist[ i ] = norm( coords1[ i ] - coords2[perm[ i ]] )


Solving the shape-matching problem using IRA
============================================

For two structures which are also rotated with respect to each other, the IRA algorithm is used
to obtain the rotation, permutation, and translation.

   >>> ## set the first atomic structure
   >>> nat1 = 10
   >>> typ1 = np.ones( nat1, dtype=int )
   >>> coords1 = np.array([[-0.49580341,  0.9708181 ,  0.37341428],
   ...                     [-1.05611656, -0.4724503 , -0.37449784],
   ...                     [-0.63509644, -0.66670776,  0.66219897],
   ...                     [-0.83642178,  0.59155936, -0.64507703],
   ...                     [ 0.59636159,  0.80558701,  0.23843962],
   ...                     [ 0.25975284,  0.71540297, -0.78971024],
   ...                     [-0.09743308, -1.03812804, -0.31233049],
   ...                     [ 0.09254502,  0.20016738,  1.03021068],
   ...                     [-0.18424967, -0.24756757, -1.07217522],
   ...                     [ 0.46705991, -0.73516435,  0.56288325]])
   >>> ##
   >>> ## set the second atomic structure
   >>> nat2 = 10
   >>> typ2 = np.ones(nat2, dtype=int)
   >>> coords2 = np.array([[-1.10284703,  0.14412375,  0.19443024],
   ...                     [ 0.66659232,  0.55627796, -0.56721304],
   ...                     [ 0.48071837,  0.23696574,  1.09688377],
   ...                     [ 1.08098955, -0.07699871,  0.21481947],
   ...                     [-0.66132935, -0.27573102, -0.73025453],
   ...                     [ 0.39018548, -0.81148351,  0.65078612],
   ...                     [-0.57686949, -0.8993001 ,  0.15398734],
   ...                     [-0.42460153,  0.78820488, -0.54634801],
   ...                     [ 0.27879878,  1.07299866,  0.34477351],
   ...                     [-0.52245748, -0.27294984,  1.07110467]])
   >>> ##
   >>> ## find the shape-matching
   >>> kmax_factor = 1.8
   >>> r, t, p, hd = ira.match( nat1, typ1, coords1, nat2, typ2, coords2, kmax_factor )
   >>> ## `r` contains the rotation matrix, `t` the translation vector,
   >>> ## `p` the permutation, and ``hd`` the hasudorff distance.
   >>> hd
   0.00751228170905401

.. note::
   The ``kmax_factor`` is a multiplicative factor that needs to be larger than 1.0. Larger value increases the
   search space of the rotations, but slows down the algorithm. Default value of 1.8 seems to be quite ok.



(under construction)

ideas:

 - comparing nonequal strucs (using candidates)
 - determine kmax value
 - using some thr.



