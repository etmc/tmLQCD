include(FindPackageHandleStandardArgs)

find_library(
  TMLQCD_LEMON_LIBRARIES
  NAMES lemon
  PATH_SUFFIXES "lib" "lib64")

find_path(
  TMLQCD_LEMON_INCLUDE_DIRS
  NAMES lemon.h
  PATH_SUFFIXES "include" "include/${_pacakge_name}" "${_package_name}")

find_package_handle_standard_args(Lemon DEFAULT_MSG TMLQCD_LEMON_LIBRARIES
                                  TMLQCD_LEMON_INCLUDE_DIRS)

if(NOT TARGET tmlqcd::lemon)
  add_library(tmlqcd::lemon INTERFACE IMPORTED)
  set_target_properties(tmlqcd::lemon PROPERTIES INTERFACE_LINK_LIBRARIES
                                                 "${TMLQCD_LEMON_LIBRARIES}")
  set_target_properties(tmlqcd::lemon PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 "${TMLQCD_LEMON_INCLUDE_DIRS}")
endif()

set(TMLQCD_LEMON_FOUND ON)
mark_as_advanced(TMLQCD_LEMON_LIBRARIES TMLQCD_LEMON_INCLUDE_DIRS)
