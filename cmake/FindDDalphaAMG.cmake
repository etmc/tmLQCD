include(FindPackageHandleStandardArgs)

find_library(
  TM_DDALPHAAMG_LIBRARIES
  NAMES DDalphaAMG DDalphaAMG_devel
  PATH_SUFFIXES "lib" "lib64")

find_path(
  TM_DDALPHAAMG_INCLUDE_DIRS
  NAMES DDalphaAMG.h
  PATH_SUFFIXES "include")

find_package_handle_standard_args(
  DDalphaAMG DEFAULT_MSG TM_DDALPHAAMG_LIBRARIES TM_DDALPHAAMG_INCLUDE_DIRS)

if(TM_DDALPHAAMG_LIBRARIES
   AND TM_DDALPHAAMG_INCLUDE_DIRS
   AND NOT TARGET tmlqcd::DDalphaAMG)
  message("INCLUDE: ${TM_DDALPHAAMG_INCLUDE_DIRS}")
  add_library(tmlqcd::DDalphaAMG INTERFACE IMPORTED)
  set_property(TARGET tmlqcd::DDalphaAMG PROPERTY INTERFACE_LINK_LIBRARIES
                                                  "${TM_DDALPHAAMG_LIBRARIES}")
  set_property(
    TARGET tmlqcd::DDalphaAMG PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                       "${TM_DDALPHAAMG_INCLUDE_DIRS}")
endif()

mark_as_advanced(TM_DDALPHAAMG_LIBRARIES TM_DDALPHAAMG_INCLUDE_DIRS)
