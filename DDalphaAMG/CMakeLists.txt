# there is a lot of custom directories to circonvent the deletion of the
# CMakeLists.txt contained in the DDalphaAMG directory. CMake will clone the
# source code and build it with the default options used in the ci/cd. More
# options are available in the main CMakeLists.txt.

set(DDalphaAMG_SRC_DIR ${CMAKE_SOURCE_DIR}/DDalphaAMG/deps)

FetchContent_Declare(
  DDalphaAMG
  GIT_REPOSITORY https://github.com/etmc/DDalphaAMG.git
  SOURCE_DIR ${DDalphaAMG_SRC_DIR})

FetchContent_MakeAvailable(DDalphaAMG)

list(
  APPEND
  DDalphaAMG_SRC_GENERIC
  interpolation_generic.c
  gathering_generic.c
  sse_interpolation_generic.c
  coarse_oddeven_generic.c
  operator_generic.c
  oddeven_generic.c
  linalg_generic.c
  init_generic.c
  vcycle_generic.c
  dirac_generic.c
  coarse_operator_generic.c
  coarsening_generic.c
  schwarz_generic.c
  ghost_generic.c
  vectorization_dirac_generic.c
  linsolve_generic.c
  sse_coarse_operator_generic.c
  data_generic.c
  setup_generic.c
  sse_linalg_generic.c)

list(
  APPEND
  DDalphaAMG_HEADER_GENERIC
  interpolation_generic.h
  gathering_generic.h
  sse_interpolation_generic.h
  coarse_oddeven_generic.h
  operator_generic.h
  oddeven_generic.h
  linalg_generic.h
  init_generic.h
  vcycle_generic.h
  dirac_generic.h
  coarse_operator_generic.h
  coarsening_generic.h
  schwarz_generic.h
  ghost_generic.h
  vectorization_dirac_generic.h
  linsolve_generic.h
  sse_coarse_operator_generic.h
  data_generic.h
  setup_generic.h
  sse_linalg_generic.h
  main_pre_def_generic.h
  main_post_def_generic.h)

list(
  APPEND
  DDalphaAMG_SRC_GENERAL
  ${DDalphaAMG_SRC_DIR}/src/preconditioner.c
  ${DDalphaAMG_SRC_DIR}/src/threading.c
  ${DDalphaAMG_SRC_DIR}/src/main.c
  ${DDalphaAMG_SRC_DIR}/src/sse_dirac.c
  ${DDalphaAMG_SRC_DIR}/src/var_table.c
  ${DDalphaAMG_SRC_DIR}/src/data_layout.c
  ${DDalphaAMG_SRC_DIR}/src/linsolve.c
  ${DDalphaAMG_SRC_DIR}/src/ghost.c
  ${DDalphaAMG_SRC_DIR}/src/top_level.c
  ${DDalphaAMG_SRC_DIR}/src/dirac.c
  ${DDalphaAMG_SRC_DIR}/src/linalg.c
  ${DDalphaAMG_SRC_DIR}/src/init.c
  ${DDalphaAMG_SRC_DIR}/src/DDalphaAMG_interface.c
  ${DDalphaAMG_SRC_DIR}/src/lime_io.c
  ${DDalphaAMG_SRC_DIR}/src/sse_linalg.c
  ${DDalphaAMG_SRC_DIR}/src/solver_analysis.c
  ${DDalphaAMG_SRC_DIR}/src/io.c)

message(STATUS "${DDalphaAMG_SRC_GENERAL}")
foreach(f IN LISTS DDalphaAMG_SRC_GENERIC)
  string(REPLACE "_generic" "_float" f_float "${f}")
  message(STATUS "${f_float}")
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}.sed-done"
    COMMAND
      sed -f "${DDalphaAMG_SRC_DIR}/float.sed" "${DDalphaAMG_SRC_DIR}/src/${f}"
      > "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}"
    COMMAND ${CMAKE_COMMAND} -E touch
            "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}.sed-done"
    DEPENDS "${DDalphaAMG_SRC_DIR}/src/${f}" "${DDalphaAMG_SRC_DIR}/float.sed"
    VERBATIM)
  list(APPEND SED_MARKERS "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}.sed-done")
  list(APPEND DDalphaAMG_SRC_SINGLE_DOUBLE ${f_float})

  string(REPLACE "_generic" "_double" f_double "${f}")
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}.sed-done"
    COMMAND
      sed -f "${DDalphaAMG_SRC_DIR}/double.sed" "${DDalphaAMG_SRC_DIR}/src/${f}"
      > "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}"
    COMMAND ${CMAKE_COMMAND} -E touch
            "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}.sed-done"
    DEPENDS "${DDalphaAMG_SRC_DIR}/src/${f}" "${DDalphaAMG_SRC_DIR}/double.sed"
    VERBATIM)
  list(APPEND SED_MARKERS "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}.sed-done")
  list(APPEND DDalphaAMG_SRC_SINGLE_DOUBLE ${f_double})
endforeach()

# now parse the header
foreach(f IN LISTS DDalphaAMG_HEADER_GENERIC)
  string(REPLACE "_generic" "_float" f_float "${f}")
  message(STATUS "${f_float}")
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}.sed-done"
    COMMAND
      sed -f "${DDalphaAMG_SRC_DIR}/float.sed" "${DDalphaAMG_SRC_DIR}/src/${f}"
      > "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}"
    COMMAND ${CMAKE_COMMAND} -E touch
            "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}.sed-done"
    DEPENDS "${DDalphaAMG_SRC_DIR}/src/${f}" "${DDalphaAMG_SRC_DIR}/float.sed"
    VERBATIM)
  list(APPEND SED_MARKERS "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_float}.sed-done")
  list(APPEND DDalphaAMG_HEADER_SINGLE_DOUBLE ${f_float})

  string(REPLACE "_generic" "_double" f_double "${f}")
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}.sed-done"
    COMMAND
      sed -f "${DDalphaAMG_SRC_DIR}/double.sed" "${DDalphaAMG_SRC_DIR}/src/${f}"
      > "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}"
    COMMAND ${CMAKE_COMMAND} -E touch
            "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}.sed-done"
    DEPENDS "${DDalphaAMG_SRC_DIR}/src/${f}" "${DDalphaAMG_SRC_DIR}/double.sed"
    VERBATIM)

  list(APPEND SED_MARKERS "${CMAKE_BINARY_DIR}/DDalphaAMG/${f_double}.sed-done")
  list(APPEND DDalphaAMG_HEADER_SINGLE_DOUBLE ${f_double})
endforeach()

foreach(outfile IN LISTS DDalphaAMG_SRC_SINGLE_DOUBLE
                         DDalphaAMG_HEADER_SINGLE_DOUBLE)
  set_source_files_properties("${CMAKE_BINARY_DIR}/DDalphaAMG/${outfile}"
                              PROPERTIES GENERATED TRUE)
endforeach()

# message(STATUS "${DDalphaAMG_SRC_SINGLE_DOUBLE}")

add_custom_target(run_sed ALL DEPENDS ${SED_MARKERS})

add_library(DDalphaAMG ${DDalphaAMG_SRC_GENERAL}
                       ${DDalphaAMG_SRC_SINGLE_DOUBLE})

target_compile_options(DDalphaAMG
                       PRIVATE "$<$<COMPILE_LANG_AND_ID:C,GNU>:-mavx2;-mfma>")

add_dependencies(DDalphaAMG run_sed)

target_link_libraries(
  DDalphaAMG
  PUBLIC MPI::MPI_C $<$<BOOL:${DDalphaAMG_ENABLE_HDF5}>:hdf5:hdf5>
         $<$<BOOL:${DDalphaAMG_ENABLE_OMP}>:OpenMP::OpenMP_C> tmlqcd::clime)

target_include_directories(
  DDalphaAMG
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${DDalphaAMG_SRC_DIR}/src>
         $<BUILD_INTERFACE:${DDalphaAMG_SRC_DIR}/include>
         $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/DDalphaAMG>)

target_compile_definitions(
  DDalphaAMG
  PUBLIC
    $<$<BOOL:${DDalphaAMG_ENABLE_PARAMOUNT_OUTPUT}>:PARAMOUNTOUTPUT>
    $<$<BOOL:${DDalphaAMG_ENABLE_FGMRES_RESTEST}>:FGMRES_RESTEST>
    $<$<BOOL:${DDalphaAMG_ENABLE_PROFILING}>:PROFILING>
    $<$<BOOL:${DDalphaAMG_ENABLE_SINGLE_ALLREDUCE_ARNOLDI}>:SINGLE_ALLREDUCE_ARNOLDI>
    $<$<BOOL:${DDalphaAMG_ENABLE_COARSE_RES}>:COARSE_RES>
    $<$<BOOL:${DDalphaAMG_ENABLE_SCHWARZ_RES}>:SCHWARZ_RES>
    $<$<BOOL:${DDalphaAMG_ENABLE_OMP}>:OPENMP>
    $<$<BOOL:${DDalphaAMG_ENABLE_TRACK_RES}>:TRACK_RES>
    $<$<BOOL:${DDalphaAMG_ENABLE_TESTVECTOR_ANALYSIS}>:TESTVECTOR_ANALYSIS>
    $<$<BOOL:${DDalphaAMG_ENABLE_HDF5}>:HAVE_HDF5>
    $<$<CONFIG:Debug>:DEBUG>
    SSE)

# add_library(tmlqcd::DDalphaAMG alias DDalphaAMG)
