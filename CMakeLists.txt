cmake_minimum_required(VERSION 3.5)

## ----------- local function
## find (sub-)string in file; return whole line with string
function(get_string in_string in_file out_string)
  if(NOT EXISTS "${in_file}")
    message(
      FATAL_ERROR
      "IRA :: file ${in_file} does not exist."
      )
  endif()

  # # read file into array of strings
  file(STRINGS ${in_file} this_file)

  set(found "")
  # parse file strings
  while(this_file)
    # pop an array element into LINE
    list(POP_FRONT this_file LINE)
    # this is a regex check if `string` matches `LINE`
    if(${LINE} MATCHES ${in_string})
      set(found ${LINE})
      break()
    endif()
  endwhile()

  set(${out_string} ${found} PARENT_SCOPE)
endfunction(get_string)
## ------ end local


# get ira verion string from file
get_string("string =" ${CMAKE_SOURCE_DIR}/src/version.f90 vstr)
# match string between double quotes
string(REGEX MATCH "^.+\"(.+)\"" _match "${vstr}")
set(ira_version ${CMAKE_MATCH_1})


project(IRA
	VERSION ${ira_version}
	DESCRIPTION "Algorithms for Iterative Rotations and Assignments (IRA); Symmetry Operations FInder (SOFI); and Constrained Shortest Distance Assignment (CShDA)."
  LANGUAGES Fortran)
enable_language(Fortran)

# find lapack/blas
find_package(LAPACK REQUIRED)
## TODO: if(NOT FOUND LAPACK) compile and link local

# add directory
add_subdirectory(src)
