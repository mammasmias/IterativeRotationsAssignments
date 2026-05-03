enable_testing()


# basic program with f90 caller
add_executable( f90main main.f90 )
set_target_properties( f90main PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE )
target_compile_options( f90main PRIVATE
  -ffree-line-length-0)
target_link_libraries( f90main PRIVATE ira )
add_test( NAME f90main COMMAND f90main )

# basic program with C caller
add_executable( cmain main.c )
target_link_libraries( cmain PRIVATE ira )
add_test( NAME cmain COMMAND cmain )


enable_language(CXX)
add_executable( cppmain main.cpp )
target_link_libraries( cppmain PRIVATE ira )
add_test( NAME cppmain COMMAND cppmain )
