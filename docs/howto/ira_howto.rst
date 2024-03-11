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

The corresponding algorithm class (IRA or SOFI) has to be initialised by:

>>> ira = ira_mod.IRA()

or

>>> sofi = ira_mod.SOFI()

now, the functions in either class are availably by typing ``ira.<funtion_name>`` or ``sofi.<function_name>``.
Quick help can be accessed by ``help( ira )``, ``help( ira.<function_name> )`` or the same for ``sofi``.


.. warning::
   If the ``ira_mod`` module cannot be found at ``import``, then make sure there is a path to ``/IRA_library/interface``
   in the environment variable ``PYTHONPATH``.
   
   .. code-block:: bash
   
      echo $PYTHONPATH
   
   If not, add it by:
   
   .. code-block:: bash
   
      export PYTHONPATH=/your/path/to/IRA_library/interface:$PYTHONPATH


Direct comparison of structures using CShDA
===========================================

ideas:

 - comparing equal strucs
 - comparing nonequal strucs (using candidates)
 - determine kmax value
 - using some thr.



