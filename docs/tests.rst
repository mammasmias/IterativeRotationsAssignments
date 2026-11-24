.. _tests:

Tests and benchmarks
====================

The tests are built with ``cmake``, by configuring with ``-DIRA_BUILD_TESTS=ON``, and are run by ``ctest``:

.. code-block:: bash

   cmake -B build -DIRA_BUILD_TESTS=ON
   cmake --build build
   ctest --test-dir build


Running benchmarks
------------------

The benchmark tests use the datasets from the original IRA publication, which are not part of this repository.
Download the IRA benchmark data from `zenodo <https://zenodo.org/doi/10.5281/zenodo.10568513>`_,
and unpack it into a directory, e.g. ``bm_data``.
Each dataset in ``bm_data/data/`` is its own archive, unpack those too:

.. code-block:: bash

   mkdir -p bm_data
   tar -xzf benchmark_test.tar.gz -C bm_data
   for f in bm_data/data/*.tar.gz; do tar -xzf "$f" -C bm_data/data; done

Configure with tests and benchmarks, pointing ``-DIRA_BM_DATA_DIR`` to the ``bm_data/data`` directory
(it should contain ``<dataset>/xyz/*.xyz``), then build:

.. code-block:: bash

   cmake -B build -DIRA_BUILD_TESTS=ON -DIRA_BUILD_BENCHMARK=ON -DIRA_BM_DATA_DIR=$PWD/bm_data/data
   cmake --build build

Each dataset is a separate test, with the name ``ira_bm_<dataset>``, and the label ``benchmark``.
Run all of them, or a single one:

.. code-block:: bash

   ctest --test-dir build -L benchmark
   ctest --test-dir build -R ira_bm_Au26 -V

The available options are:

* ``-DIRA_BUILD_BENCHMARK``: build the benchmark program and tests, requires ``-DIRA_BUILD_TESTS=ON`` (default OFF);
* ``-DIRA_BM_DATA_DIR``: path to the benchmark datasets (default ``IRA/benchmark_test/ira_bm/data``);
* ``-DIRA_BM_INCLUDE_LARGE``: include the ``lj_clusters`` dataset, which contains 859 structures (default OFF).

.. admonition:: running only the basic tests
   :class: tip

   When the benchmarks are configured, they are also run by a plain ``ctest``. To exclude them, do:

   .. code-block:: bash

      ctest --test-dir build -LE benchmark


The ``ira_bm`` program
----------------------

The benchmarking program is called ``ira_bm``, and is located in ``build/tests/IRA_bm/``.
Standalone, it takes ``.xyz`` files as arguments:

.. code-block:: bash

   build/tests/IRA_bm/ira_bm file1.xyz file2.xyz

The flags are:

* ``--list <file>``: read the ``.xyz`` filenames from ``<file>``, one path per line. The lists for each dataset are written
  by ``cmake`` to ``build/tests/IRA_bm/<dataset>.list``;
* ``--log-file``: print the time spent per structure;
* ``--log-each``: print details of every match.

For example:

.. code-block:: bash

   build/tests/IRA_bm/ira_bm --log-file --list build/tests/IRA_bm/Au26.list
