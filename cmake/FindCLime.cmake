include(FindPackageHandleStandardArgs)

find_library(
  TM_CLIME_LIBRARIES
  NAMES lime
  PATH_SUFFIXES "lib" "lib64")

find_path(
  TM_CLIME_INCLUDE_DIRS
  NAMES lime.h
  PATH_SUFFIXES "include" "include/${_pacakge_name}" "${_package_name}")

find_package_handle_standard_args(CLime DEFAULT_MSG TM_CLIME_LIBRARIES
                                  TM_CLIME_INCLUDE_DIRS)

if(NOT TARGET tmlqcd::clime)
  add_library(tmlqcd::clime INTERFACE IMPORTED)
  set_target_properties(tmlqcd::clime PROPERTIES INTERFACE_LINK_LIBRARIES
                                                 "${TM_CLIME_LIBRARIES}")
  set_target_properties(tmlqcd::clime PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 "${TM_CLIME_INCLUDE_DIRS}")
endif()

set(TM_CLIME_FOUND ON)
mark_as_advanced(TM_CLIME_FOUND TM_CLIME_LIBRARIES
                 TM_CLIME_INCLUDE_DIRS)
