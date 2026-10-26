# Description

This is the repository of shape matching algorithm
Iterative Rotations and Assignments (IRA), and the Symmetry Operations FInder (SOFI) algorithm.

Online documentation and tutorials: [LINK](https://mammasmias.github.io/IterativeRotationsAssignments/).
GitHub repository: [LINK](https://github.com/mammasmias/IterativeRotationsAssignments).

IRA is described in the publication [[1]](#1).
It is also the main subject of the dissertation [[2]](#2), where a workflow
inserting IRA into an off-lattice kMC algorithm is developed.

SOFI is described in the publication [[3]](#3).

# Directory contents

 `/src`:
   Contains the source code of the algorithms, and the C-bound API library.

 `/examples`:
   Contains example programs which use different functionalities of the IRA/CShDA/SOFI algorithms.

 `/interface`:
   Contains the C-headers, and Python module to interface the API.

 `/benchmark_test`:
   Contains data and other software used for benchmark tests done in [[1]](#1). See
   also `/benchmark_test/README`. NOTE: the contents have been moved to [zenodo](https://zenodo.org/doi/10.5281/zenodo.10568513).


# Compile and run

For the impatient:

    python -m pip install .

Python module:

    >>> import ira_mod
    >>> ira_mod.version
    >>> ira = ira_mod.IRA()
    >>> sofi = ira_mod.SOFI()

For other build tools and more details, see the [online documentation](https://mammasmias.github.io/IterativeRotationsAssignments/compilation.html).

The shared library hides the reference BLAS that the bundled `lap.f` would otherwise export. A caller that used those symbols (`dgemm_`, `dgemv_`, `dgesvd_`, `xerbla_`) from `libira` must link a BLAS itself. Opt out with `-DIRA_HIDE_BUNDLED_BLAS=OFF` or `IRA_HIDE_BUNDLED_BLAS=0`. GNU ld applies `src/libira.map`. A linker that rejects `-Wl,--version-script` warns and hides the same symbols with `-Wl,--exclude-libs` on the bundled archive. If that flag is rejected too, the shared build stops rather than exporting `dgemm_`.

# Conda build

Using the conda-forge channel:

    conda install ira


# Terms and conditions
The software in this repository is subject to the license(s) provided in the `LICENSE.txt` file.


## References

<a id="1">[1]</a>
Gunde M., Salles N., Hemeryck A., Martin Samos L.
*IRA: A shape matching approach for recognition and comparison of generic atomic patterns*,
Journal of Chemical Information and Modeling (2021), DOI:
[https://doi.org/10.1021/acs.jcim.1c00567](https://doi.org/10.1021/acs.jcim.1c00567),
HAL: [hal-03406717](https://laas.hal.science/hal-03406717v1), arXiv:
[2111.00939](https://export.arxiv.org/abs/2111.00939)

<a id="2">[2]</a>
Gunde M.: *Development of IRA: a shape matching algorithm, its implementation
and utility in a general off-lattice kMC kernel*, PhD dissertation, Université Toulouse III - Paul Sabatier,
November 2021,
[link](https://theses.hal.science/tel-03635139v2).

<a id="3">[3]</a>
Gunde M., Salles N., Grisanti L., Hemeryck A., Martin Samos L.
*SOFI: Finding point group symmetries in atomic clusters as finding the set of degenerate solutions in a shape-matching problem*,
Journal of Chemical Physics (2024), DOI:
[https://doi.org/10.1063/5.0215689](https://doi.org/10.1063/5.0215689),
arXiv: [2408.06131](https://arxiv.org/abs/2408.06131).
