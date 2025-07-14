.. _compilation:

Compilation
===========

Traditional ``make``
--------------------

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


.. admonition:: lapack library
   :class: tip

   The lapack library is needed for compilation.
   Default flag is ``LIBLAPACK=-llapack``, however if you wish to change that, i.e. with openblas, you can specify it from the command as:

   .. code-block:: bash

      LIBLAPACK=-lopenblas make all

   If the lapack library is not available on your system, you can leave the variable undefined (this will compile a local version of the needed lapack routines, which is however not optimal):

   .. code-block:: bash

      LIBLAPACK='' make all


``conda`` compatible build
--------------------------

For local machines, it is possible to use ``pixi`` [pixi-install]_ to get a working version of the
Python bindings in a fairly automated manner.

.. code-block:: bash

   curl -fsSL https://pixi.sh/install.sh | bash
   # build the library with openblas
   pixi run build_lib
   pixi shell # sets the environment variable
   cd examples/IRA
   python python_program.py

.. [pixi-install] Installation instructions here: `<https://pixi.sh/latest/>`_



Using ``cmake``
---------------

To install with ``cmake``, it is assumed you have the ``lapack`` or ``blas`` library installed on your system.
It will install the library into ``IRA/lib`` directory, and the module file into ``IRA/include``, regardless of the ``CMAKE_INSTALL_PREFIX`` variable.

.. code-block:: bash

   cmake -B build
   cmake --build build
   cmake --install build



Linking a program to libira
===========================

A program compiled with ``gcc`` or ``gfortran`` can easily link the IRA library, as-is, by linking either the shared
library ``libira.so``, or the static version ``libira.a``. They are both located in the ``lib/`` directory after
compilation.

Example for fortran program:

.. code-block:: bash

   gfortran -o caller_program.x caller_program.f90 -L/your/path/to/IRA/lib/ -lira -Wl,-rpath,/your/path/to/IRA/lib

The base-level implementations are not placed in modules, therefore all routines are in principle acessible to the
caller. Care must be taken to ensure the correct type, kind, shape, etc. of the arguments, i.e. interface matching
needs to be checked manually.
The default precision is equivalent to ``c_int`` for integers, and ``c_double`` for reals, they are defined in ``IRA/src/ira_precision.f90`` module.

The C-headers are located in the ``IRA/interface`` directory, and can be included in compilation by ``-I/your/path/to/IRA/interface``.

When linking the static library ``libira.a`` to a C-program, you need to add the math (``-lm``), and fortran (``-lgfortran``, or equivalent) to the compilation:

.. code-block:: bash

   gcc -I/your/path/IRA/interface -o c_prog.x c_prog.c -L/your/path/to/IRA/src -lira -Wl,-rpath,/your/path/to/IRA/src -lm -lgfortran


