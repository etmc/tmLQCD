include(FindPackageHandleStandardArgs)

find_library(
  TM_DDALPHAAMG_LIBRARIES
  NAMES DDalphaAMG DDalphaAMG_devel
  PATH_SUFFIXES "lib" "lib64")

find_path(
  TM_DDALPHAAMG_INCLUDE_DIRS
  NAMES DDalphaAMG.h
  PATH_SUFFIXES "include" "include/${_pacakge_name}" "${_package_name}")

find_package_handle_standard_args(
  DDAlphaAMG DEFAULT_MSG TMLQCD_DDALPHAAMG_LIBRARIES
  TMLQCD_DDALPHAAMG_INCLUDE_DIRS)

if(NOT TARGET tmlqcd::DDalphaAMG)
  add_library(tmlqcd::DDalphaAMG INTERFACE IMPORTED)
  set_target_properties(
    tmlqcd::DDalphaAMG PROPERTIES INTERFACE_LINK_LIBRARIES
                                  "${TMLQCD_DDALPHAAMG_LIBRARIES}")
  set_target_properties(
    tmlqcd::DDalphaAMG PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                  "${TMLQCD_DDALPHAAMG_INCLUDE_DIRS}")
endif()

set(TMLQCD_DDALPHAAMG_FOUND ON)
mark_as_advanced(TMLQCD_DDALPHAAMG_FOUND TMLQCD_DDALPHAAMG_LIBRARIES
                 TMLQCD_DDALPHAAMG_INCLUDE_DIRS)
