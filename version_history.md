Version history:

 * Version 1.0.0: November 2021

 the IRA code is released without the data and scripts of the benchmark tests
 from the original publication.


 * Version 1.5.0: May 2022

 changes with respect to previous version:
  - added a threshold for early exit of the cshda computation;
  - the threshold from the latter point is updated within the main IRA loop in a self-consistent manner, which significantly improves the speed of IRA;
  - unification of eq and noneq routines into ira_unify;
  - added set_candidate routines for modifiable control of candidate central atoms;
  - added example calling program from f90;
  - added interface to python in ira_interf.pyf, which compiles into a python module using f2py3;
  - added example calling program from python.


 * Version 2.0.0: April 2024 (20240425)

 Major changes with respect to previous version:
  - the SOFI algorithm for finding point group symmetries is added;
  - all routines have a BIND(C) equivalent, the library is fully interoperable with C;
  - replaced f2py with ctypes for python interfacing;
  - added routine for obtaining the current version, get_version();
  - online documentation with some tutorials is added.


 * Version 2.1.0: July 2024 (20240718)

 Minor changes:
  - add functionality to deal with linear structures in SOFI;
  - add some example inputs for SOFI;
  - some changes in API calls;
  - slight speedup in CShDA achieved through more aggressive compiler flags;


 * Version 2.2.0: September 2025 (20250911)

 Minor changes:
  - add `ira_precision`;
  - add more build systems (pip, pixi, cmake, fpm);
  - add local version of lapack routines if needed;
  - solve a rare bug in cshda when different number of atoms;
