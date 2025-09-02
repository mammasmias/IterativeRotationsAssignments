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


# Terms and conditions
The software in this repository is subject to the license(s) provided in the `LICENSE.txt` file.


## References

<a id="1">[1]</a>
Gunde M., Salles N., Hemeryck A., Martin Samos L.
*IRA: A shape matching approach for recognition and comparison of generic atomic patterns*,
Journal of Chemical Information and Modeling (2021), DOI:
[https://doi.org/10.1021/acs.jcim.1c00567](https://doi.org/10.1021/acs.jcim.1c00567),
HAL: [hal-03406717](https://hal.laas.fr/hal-03406717), arXiv:
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
