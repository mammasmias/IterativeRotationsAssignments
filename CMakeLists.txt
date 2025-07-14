cmake_minimum_required(VERSION 3.5)

project(IRA 
	VERSION 2.1.0
	DESCRIPTION "Algorithms for Iterative Rotations and Assignments (IRA); Symmetry Operations FInder (SOFI); and Constrained Shortest Distance Assignment (CShDA)."
        LANGUAGES Fortran)
enable_language(Fortran)

# set obj dir for module files (slash at end is important for install path it seems)
# CMAKE_BINARY_DIR is the dir specified at `cmake -B <dir>`
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/Obj/)

# compiler flags
add_compile_options(-cpp -ffree-line-length-none -march=native -funroll-loops -ffast-math -Ofast)

# find lapack/blas
find_package(LAPACK REQUIRED)

# add directory
add_subdirectory(src)



## targets: these get executed at `cmake --install <dir>`
## the paths are by default into: CMAKE_INSTALL_PREFIX/<DESTINATION>
## where the <DESTINATION> is specified below.
## However we need the lib/ dir to be fixed due to python. So force the destination:
## CMAKE_CURRENT_SOURCE_DIR is where the root CMakeLists is located.

# install target: library libira.so
install(TARGETS ira 
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# install target: module file
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include
	FILES_MATCHING PATTERN "ira_precision.mod" )


