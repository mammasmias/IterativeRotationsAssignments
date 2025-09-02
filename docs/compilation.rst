.. _compilation:

Compilation
===========

The compilation will produce directories ``IRA/lib`` and ``IRA/include``, which contain the library and module files respectively. Other directories may be creted, depending on the flavour of build tool used.

Python module via ``pip``
-------------------------

If you just want the ``ira_mod`` python module, it can be done by running in ``IRA/`` root directory (this relies on ``cmake``):

.. code-block:: bash

   python -m pip install .

Then you can directly use it:

.. code-block:: python

   >>> import ira_mod
   >>> print( ira_mod.version )
   >>> ira = ira_mod.IRA()
   >>> sofi = ira_mod.SOFI()


Other build tools
-----------------

Other build tools are available, use one of the following:

.. tab-set::

   .. tab-item:: Traditional ``make``

      To compile the IRA library, you need the ``lapack`` library.
      On a standard linux machine, it should suffice to type:

      .. code-block:: bash

         cd src/
         make all


      This will generate the static and shared libraries: ``libira.a``, and ``libira.so``.
      To compile only one of the libraries, type ``make lib`` or ``make shlib``, respectively.

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


      .. admonition:: python module ``ira_mod``
         :class: tip

         To use the python module, you will need to set the ``PYTHONPATH`` variable:

         .. code-block::

            export PYTHONPATH=/path/to/IRA/interface:$PYTHONPATH


   .. tab-item:: ``conda`` compatible build

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


   .. tab-item:: Using ``cmake``

      To install with ``cmake``, it is assumed you have the ``lapack`` or ``blas`` library installed on your system.

      .. code-block:: bash

         cmake -B builddir
         cmake --build builddir

      .. admonition:: python module ``ira_mod``
         :class: tip

         To use the python module, you will need to set the ``PYTHONPATH`` variable:

         .. code-block::

            export PYTHONPATH=/path/to/IRA/interface:$PYTHONPATH



   .. tab-item:: Using ``fpm``

      Required minimum ``fpm`` version 0.12.0.

      .. code-block:: bash

         fpm build --flag "-fPIC -fcheck=bounds -ffree-line-length-none -Ofast -march=native -ffast-math -funroll-loops"
         fpm install --prefix .

      .. admonition:: python module ``ira_mod``
         :class: tip

         To use the python module, you will need to set the ``PYTHONPATH`` variable:

         .. code-block::

            export PYTHONPATH=/path/to/IRA/interface:$PYTHONPATH



Linking a program to libira
===========================

A program compiled with ``gcc`` or ``gfortran`` can easily link the IRA library, as-is, by linking either the shared
library ``libira.so``, or the static version ``libira.a``. They are both located in the ``lib/`` directory after
compilation. The module files are located in ``include/``.

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


