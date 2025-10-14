include(FindPackageHandleStandardArgs)

find_library(
  TMLQCD_CLIME_LIBRARIES
  NAMES lime
  PATH_SUFFIXES "lib" "lib64")

find_path(
  TMLQCD_CLIME_INCLUDE_DIRS
  NAMES lime.h
  PATH_SUFFIXES "include" "include/${_pacakge_name}" "${_package_name}")

message("${TMLQCD_CLIME_INCLUDE_DIRS}")
find_package_handle_standard_args(CLime DEFAULT_MSG TMLQCD_CLIME_LIBRARIES
                                  TMLQCD_CLIME_INCLUDE_DIRS)

if(NOT TARGET tmlqcd::clime)
  add_library(tmlqcd::clime INTERFACE IMPORTED)
  set_target_properties(tmlqcd::clime PROPERTIES INTERFACE_LINK_LIBRARIES
                                                 "${TMLQCD_CLIME_LIBRARIES}")
  set_target_properties(tmlqcd::clime PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 "${TMLQCD_CLIME_INCLUDE_DIRS}")
endif()

set(TMLQCD_CLIME_FOUND ON)
mark_as_advanced(TMLQCD_CLIME_FOUND TMLQCD_CLIME_LIBRARIES
                 TMLQCD_CLIME_INCLUDE_DIRS)
