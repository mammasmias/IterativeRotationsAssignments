.. _ref_api:

###############
The C-bound API
###############

The API is defined in files ``library_ira.f90`` and ``library_sofi.f90``. It can be used
to call IRA/SOFI routines directly from C. Calling from other languages needs an interface, see
for example :ref:`py_interf`.

.. note::

    This file defines the wrappers to IRA routines, which are part of the
    shared library libira.so.
    The routines here are defined with "bind(C)" and "use iso_c_binding", and
    they effectively behave as C-code. **The arrays passed to these routines are
    assumed to be already allocated by the caller**, and are assumed to have
    C-shape. Take care of transposed matrices, starting indices (1 or 0), etc.
    Likewise, the **output arrays are assumed to be already allocated by the caller,
    and it's up to the caller to assure the correct shape and indices**.
    The routines here receive C-pointers to the memory where output should be
    written. Therefore, **the output data appears as "intent(in)"**.

.. note::

   There is no automatic checking of the type/kind/shape of arguments in the API.
   The **programmer is trusted** that her/his calls to the API functions are correct.

.. contents:: Contents
   :local:
   :depth: 1


CShDA & IRA
===========

.. doxygenfile:: library_ira.f90
   :project: lib_IRA


SOFI
====

.. note::

   The size of arguments in many of the SOFI functions are pre-defined with the ``nmax`` variable, which
   is defined in sofi_tools.f90, and is by default ``nmax=400``.
   The actual output is written to the first ``n_mat`` elements, where ``n_mat`` is the number of matrices/symmtry operations.

.. doxygenfile:: library_sofi.f90
   :project: lib_IRA
