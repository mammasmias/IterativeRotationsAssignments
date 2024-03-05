
=========================
IRA library documentation
=========================


.. contents:: Contents
   :local:
   :depth: 2

.. toctree::
   :maxdepth: 1
   :caption: Reference

   src/refs
   src/c_wrappers
   interf/py_interf




Introduction
============

IRA
---

SOFI
----


Installation
============

To compile the IRA library, you need the ``lapack`` library.
On a standard linux machine, it should suffice to type:

.. code-block:: bash

   cd src/
   make all


This will generate the static and shared libraries: ``libira.a``, and ``libira.so``.
To compile only one of them, type ``make lib`` or ``make shlib``, respectively.

To clean the build type:

.. code-block:: bash

   make clean


For customising the default compilation, see the ``src/Makefile``.
If you need to specify a custom ``lapack`` location, change the default value ``LIBLAPACK = -llapack`` to your value.


Linking a program to libira
===========================


Examples
========

