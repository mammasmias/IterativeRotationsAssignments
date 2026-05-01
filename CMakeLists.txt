cmake_minimum_required(VERSION 3.21)

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/version.txt vv )
list(GET vv 0 sem_version)
list(GET vv 1 date_version)

project(IRA
  VERSION "${sem_version}"
	DESCRIPTION
  "Algorithms for Iterative Rotations and Assignments (IRA); Symmetry Operations FInder (SOFI); and Constrained Shortest Distance Assignment (CShDA)."
  LANGUAGES Fortran C)


message(STATUS "cmake_project_version: ${CMAKE_PROJECT_VERSION}")


# add directory
add_subdirectory(src)

# ######################################################################################
# Install rules for Python packaging (scikit-build-core)
install(TARGETS ira LIBRARY DESTINATION ira_mod/lib)
install(FILES
    interface/__init__.py
    interface/ira_mod.py
    DESTINATION ira_mod
)
# ######################################################################################

