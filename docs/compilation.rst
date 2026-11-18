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


``conda`` installation
----------------------

Using the conda-forge channel, you can download and install the ``ira_mod`` python module by:

.. code-block:: bash

   conda install ira


Other build tools
-----------------

Other build tools are available, use one of the following:

.. tab-set::


   .. tab-item:: Using ``cmake`` (recommened)

      To install with ``cmake``, execute:

      .. code-block:: bash

         cmake -B builddir
         cmake --build builddir

      This will produce directories ``lib`` and ``include`` in the builddir path specified, they will contain
      the static library ``libira.a`` and the fortran module files, respectively.
      If you wish to compile the shared library, configure the build with ``-DBUILD_SHARED_LIBS=ON``:

      .. code-block:: bash

         cmake -B builddir -DBUILD_SHARED_LIBS=ON
         cmake --build builddir

      The available options are:

      * ``-DIRA_BUILD_TESTS``: build some simple tests (default OFF);
      * ``-DIRA_ENABLE_NATIVE_OPTIMIZATION`` optimize for the build machine
        (for example ``-march=native``; default ON);
      * ``-DIRA_ENABLE_FAST_MATH`` allow non-IEEE transformations (for example
        ``-Ofast -ffast-math``; default ON). Set it to ``OFF`` for ``-O3``
        without fast-math;
      * ``-DIRA_LINALG_PROVIDER=AUTO|FLEXIBLAS|LAPACK|BUNDLED`` selects the
        linear algebra provider. ``AUTO`` (the default) first asks CMake's
        ``FindBLAS`` and ``FindLAPACK`` modules for ``BLA_VENDOR=FlexiBLAS``,
        then tries any BLAS/LAPACK, and finally uses the bundled fallback.
        The other values require the selected provider.

      FlexiBLAS is preferred because one public BLAS interface can dispatch to
      a backend selected at runtime, for example with
      ``FLEXIBLAS=OPENBLASPTHREAD`` or ``FLEXIBLAS=NETLIB``. The bundled
      ``lap_local/lap.f`` is compiled only for the fallback, with ``-O2`` and
      without fast-math. Its routines are prefixed ``ira_``; ``libira`` never
      defines the public ``dgemm_`` or other BLAS/LAPACK names.

      To install the library and module files to a particular location, do:

      .. code-block:: bash

         cmake -B builddir -DCMAKE_INSTALL_PREFIX=your/install/path
         cmake --build builddir
         cmake --install builddir

      This will copy the ``lib`` and ``include`` directories and their contents to ``your/install/path``, and create files needed for other ``CMake`` and/or ``pkg-config`` based projects to find the library.

      .. admonition:: python module ``ira_mod``
         :class: tip

         To use the python module, you will need to pass the ``-DIRA_INSTALL_PYTHON=ON``, install the
         library somewhere, and set the ``PYTHONPATH`` variable to that same path with adding ``/ira_mod``:

         .. code-block::

            cmake -B builddir -DIRA_INSTALL_PYTHON=ON -DCMAKE_INSTALL_PREFIX=your/install/path
            cmake --build builddir
            cmake --install builddir
            export PYTHONPATH=your/install/path/ira_mod



   .. tab-item:: Traditional ``make``

      On a standard linux machine, it should suffice to type:

      .. code-block:: bash

         cd src/
         make all


      This will generate the static and shared libraries: ``libira.a``, and ``libira.so``.
      ``make`` asks ``pkg-config`` for FlexiBLAS first and for LAPACK second.
      If neither is found it compiles ``lap_local/lap.f`` (routines prefixed
      ``ira_``, ``-O2``). Override ``LIBLAPACK`` to select another provider;
      use ``make LIBLAPACK=`` to force the bundled fallback.
      To compile only one of the libraries, type ``make lib`` or ``make shlib``, respectively.
      ``make test-link-blas`` checks that ``libira.so`` exports none of the unprefixed BLAS names and that a program linked with ``-lira -lopenblas`` runs.

      To clean the build type:

      .. code-block:: bash

         make clean


      .. admonition:: python module ``ira_mod``
         :class: tip

         To use the python module, you will need to set the ``PYTHONPATH`` variable:

         .. code-block::

            export PYTHONPATH=/path/to/IRA/interface:$PYTHONPATH


   .. tab-item:: Using ``pixi``

      For local machines, it is possible to use ``pixi`` [pixi-install]_ to get a working version of the
      Python bindings in a fairly automated manner.

      .. code-block:: bash

         curl -fsSL https://pixi.sh/install.sh | bash
         # build the library
         pixi run build_lib
         pixi shell # sets the environment variable
         cd examples/IRA
         python python_program.py

      The conda-forge environment supplies OpenBLAS and uses the generic LAPACK
      path (``pixi run build_fpm`` or ``pixi run build_cmake``). FlexiBLAS is
      not currently packaged on conda-forge, so the Pixi environment cannot
      enable the preferred provider without a separately installed FlexiBLAS.

      .. [pixi-install] Installation instructions here: `<https://pixi.sh/latest/>`_


   .. tab-item:: Using ``fpm``

      Required minimum ``fpm`` version 0.13.0.

      .. code-block:: bash

         fpm build --profile release
         fpm install --prefix .

      The ``release`` profile links FlexiBLAS and applies the native
      ``-Ofast`` profile. ``lapack-release`` uses generic BLAS/LAPACK instead;
      ``strict-release`` keeps ``-O3 -march=native -funroll-loops`` but omits
      fast-math; and ``portable-release`` also omits ``-march=native``. Use
      ``bundled-release`` for an optimized no-dependency build or
      ``bundled-safe-release`` to compile the entire package with ``-O2``.

      fpm 0.13 applies one compiler flag set to the complete dependency graph;
      it cannot give a path dependency separate per-file flags. Consequently
      ``bundled-release`` compiles both IRA and the prefixed fallback with
      ``-O3`` (without fast-math), while ``bundled-safe-release`` keeps the
      fallback's conservative ``-O2`` at the cost of compiling IRA itself with
      the same flags. The CMake and Make builds do not have this limitation.

      .. admonition:: python module ``ira_mod``
         :class: tip

         To use the python module, you will need to set the ``PYTHONPATH`` variable:

         .. code-block::

            export PYTHONPATH=/path/to/IRA/interface:$PYTHONPATH



Linking a program to libira
===========================

Linking a program to ``libira`` is simple, however note that:

* The base-level implementations are not placed in modules, therefore all routines are in principle acessible to the
  caller. Care must be taken to ensure the correct type, kind, shape, etc. of the arguments, i.e. interface matching
  needs to be checked manually;
* The default precision is equivalent to ``c_int`` for integers, and ``c_double`` for reals, they are defined in
  the ``IRA/src/ira_precision.f90`` module.


.. tab-set::
   .. tab-item:: Using ``cmake`` (recommended)

      To link your Fortran or C target in CMake with ``libira``, it suffices to add IRA project
      into your ``CMakeLists.txt`` project by one of the following methods:

      * To include IRA as part of the build process of your program:

        * With ``FetchContent`` as:

          .. code-block:: cmake

             # add IRA via FetchContent:
             include( FetchContent )
             FetchContent_Declare( ira_git
               GIT_REPOSITORY https://github.com/mammasmias/IterativeRotationsAssignments.git
               GIT_TAG master
               GIT_SHALLOW 1
               )
             FetchContent_MakeAvailable( ira_git )

             # link it to your target:
             target_link_libraries( your-target PRIVATE ira::ira )

        * With ``add_subdirectory()`` as submodule as:

          .. code-block:: cmake

             # i.e. IRA is a submodule located in a directory 'path/to/IRA'
             add_subdirectory( path/to/IRA )

             # link it to your target:
             target_link_libraries( your-target PRIVATE ira::ira )

      * To include an already-compiled IRA somewhere on your machine:

        * With ``include()`` to locate the exported targets file. Note, IRA exports a config file during the build step, located in the binary/cmake/ira dir, and during install step located in the install/lib/cmake/ira dir.

          .. code-block:: cmake

             # using include() with the build directory:
             include( path-to-ira-build/cmake/ira/ira-BuildConfig.cmake )

             # using include() with the install directory:
             include( path-to-ira-install/lib/cmake/ira/ira-config.cmake )

             # link to your target:
             target_link_libraries( your-target PRIVATE ira::ira )

        * Or with ``find_package()`` as:

          .. code-block:: cmake

             # using find_package() with the build directory:
             find_package( ira HINTS path-to-ira-build/cmake/ira NAMES ira-Build )

             # using find_package() with the install directory:
             find_package( ira HINTS path-to-ira-install )

             # link to your target:
             target_link_libraries( your-target PRIVATE ira::ira )


   .. tab-item:: Using traditional ``make``

      A program compiled with ``gcc`` or ``gfortran`` can easily link the IRA library, as-is, by linking either the shared
      library ``libira.so``, or the static version ``libira.a``. They are both located in the ``lib/`` directory after
      compilation. The module files are located in ``include/``.

      Example for fortran program:

      .. code-block:: bash

         gfortran -o caller_program.x caller_program.f90 -L/your/path/to/IRA/lib/ -lira -Wl,-rpath,/your/path/to/IRA/lib

      The C-headers are located in the ``IRA/interface`` directory, and can be included in compilation by ``-I/your/path/to/IRA/interface``.

      When linking the static library ``libira.a`` to a C-program, you need to add the math (``-lm``), and fortran (``-lgfortran``, or equivalent) to the compilation:

      .. code-block:: bash

         gcc -I/your/path/IRA/interface -o c_prog.x c_prog.c -L/your/path/to/IRA/src -lira -Wl,-rpath,/your/path/to/IRA/src -lm -lgfortran


