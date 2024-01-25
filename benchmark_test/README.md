=======================
NOTE:
=====

The contents of this directory have been moved to: [zenodo](https://zenodo.org/doi/10.5281/zenodo.10568513). 
The original `README.md` file is kept below, for the sake of completeness.

=======================





# Directory contents
Benchmark test of the IRA algorithm has been done in [[1]](#1). Present directory
contains the data and software, other than IRA, used in the test.

Description of directory contents:

 `/ArbAlign`:
   Contains the ArbAlign software, obtained from [[2]](#2) with git repository [https://github.com/berhane/arbalign](https://github.com/berhane/arbalign). See also README in
   the directory.

 `/convert_xyz2FO`:
   Contains the converter from xyz to format specific to fastoverlap.

 `/data`:
   Contains data used for the benchmark test. See also README
   in the directory. The data has been obtained from other published works
   [[2]](#2)[[3]](#3)[[4]](#4)[[5]](#5)[[6]](#6)

 `/fastoverlap`:
   Contains the fastoverlap software, obtained from [[7]](#7) with git repository [https://github.com/matthewghgriffiths/fastoverlap](https://github.com/matthewghgriffiths/fastoverlap). See also README
   in the directory.

 `/randomize`:
   Contains the software to randomize a structure, see README therein.

# Compile and run:

To run `/randomize`, and `/convert_xyz2FO`, you need to first compile them.
For that, see the README files in appropriate subdirectory.
To run `/fastoverlap` you need to unzip and copy a file, see `/fastoverlap/README`.

# Scripts

Description of scripts:
 `run_single.sh`:
   Runs a single instance of the randomization test. See header of the file.

 `run_many.sh`:
   Runs the randomization test over a directory of data. See header of the file.


## Terms and conditions

The ArbAlign and fastoverlap software, used for benchmark test purposes in the
`/benchmark_test` directory,
are subject to their own licenses. According to their original publications,
they are GNU-GPL v2 for ArbAlign [[2]](#2), and GNU-GPL v3 for fastoverlap
[[7]](#7).

The data used for benchmark test are originated from several sources, please refer to
the original publicaions mentioned in the README files
(`/benchmark_test/data/*/README`) of each data set for their terms and conditions.


## References

<a id="1">[1]</a> Gunde M., Salles N., Hemeryck A., Martin Samos L. *IRA: A
shape matching approach for recognition and comparison of generic atomic
patterns*, Journal of Chemical Information and Modeling (2021), DOI:
[https://doi.org/10.1021/acs.jcim.1c00567](https://doi.org/10.1021/acs.jcim.1c00567),
HAL: [hal-03406717](https://hal.laas.fr/hal-03406717), arXiv:
[2111.00939](https://export.arxiv.org/abs/2111.00939); and Gunde M.:
*Development of IRA: a shape matching algorithm, its implementation and utility
in a general off-lattice kMC kernel*, PhD dissertation, November 2021. [PDF
link](http://thesesups.ups-tlse.fr/5109/1/2021TOU30132.pdf)

<a id="2"> [2]</a> B. Temelso, J. M. Mabey, T. Kubota, N. Appiah-Padi, G. C.
Shields, J. Chem. Inf. Model. 2017 57, 1054, DOI:
[10.1021/acs.jcim.6b00546](https://doi.org/10.1021/acs.jcim.6b00546)

<a id="3">[3]</a> Shao X., Wu X., Cai W., J. Phys. Chem. A, 2010 114, 29,
DOI: [10.1021/jp906922v](https://doi.org/10.1021/jp906922v)

<a id="4">[4]</a> Liu Q., Xu C., Wu X., Cheng L., Nanoscale 2019 11, 13227, DOI:
[10.1039/C9NR02617G](http://dx.doi.org/10.1039/C9NR02617G)

<a id="5">[5]</a> Brena B., Ojamäe L., J. Phys. Chem. C 2008 112,
13516, DOI: [10.1021/jp8048179](https://doi.org/10.1021/jp8048179)

<a id="6">[6]</a> D. J. Wales, J. P. K. Doye, A. Dullweber, M. P.
Hodges, F. Y. Naumkin, F. Calvo, J. Hern&agrave;ndez-Rojas, T. F.
Middleton, The Cambridge Cluster Database, URL:[https://www-wales.ch.cam.ac.uk/CCD.html](https://www-wales.ch.cam.ac.uk/CCD.html)

<a id="7">[7]</a> Griffiths M. Niblett S. P., Wales D. J., J. Chem. Theory
Comput. 2017 13, 4914, DOI:
[10.1021/acs.jctc.7b00543](https://doi.org/10.1021/acs.jctc.7b00543)
