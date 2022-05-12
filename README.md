# Description

This is the repository of shape matching algorithm 
Iterative Rotations and Assignments (IRA), described in the publication
[[1]](#1). It is also the main subject of the dissertation [[2]](#2), where a workflow
inserting IRA into an off-lattice kMC algorithm is developed.


# Directory contents

 `/IRA`: 
   Contains the IRA software, see also `/IRA/README`.

 `/benchmark_test`:
   Contains data and other software used for benchmark tests done in [[1]](#1). See
   also `/benchmark_test/README`.


# Compile and run
To run IRA, you need to compile it. See `/IRA/README`.


# Terms and conditions
The IRA software, contained in the `/IRA` directory, is subject to the license(s)
provided in the `/IRA/LICENSE.txt` file.

The ArbAlign and fastoverlap software, used for benchmark test purposes in the
`/benchmark_test` directory,
are subject to their own licenses. According to their original publications,
they are GNU-GPL v2 for ArbAlign [[3]](#3) with the git repository at [https://github.com/berhane/arbalign](https://github.com/berhane/arbalign), and GNU-GPL v3 for fastoverlap
[[4]](#4) with the git repository at
[https://github.com/matthewghgriffiths/fastoverlap](https://github.com/matthewghgriffiths/fastoverlap).

The data used for benchmark test are originated from several sources, please refer to
the original publicaions mentioned in the README files
(`/benchmark_test/data/*/README`) of each data set for their terms and conditions.

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
and utility in a general off-lattice kMC kernel*, PhD dissertation,
November 2021.
[PDF link](http://thesesups.ups-tlse.fr/5109/1/2021TOU30132.pdf) 

<a id="3"> [3]</a> B. Temelso, J. M. Mabey, T. Kubota, N. Appiah-Padi, G. C.
Shields, J. Chem. Inf. Model. 2017 57, 1054, DOI:
[10.1021/acs.jcim.6b00546](https://doi.org/10.1021/acs.jcim.6b00546)

<a id="4">[4]</a> Griffiths M. Niblett S. P., Wales D. J., J. Chem. Theory
Comput. 2017 13, 4914, DOI:
[10.1021/acs.jctc.7b00543](https://doi.org/10.1021/acs.jctc.7b00543)
