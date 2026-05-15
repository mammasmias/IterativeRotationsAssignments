@PACKAGE_INIT@

# if the target ira::ira does not exist, include the targets file to create the targets
if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@" )
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
endif()
