list(APPEND tmlqcd_prog "benchmark;deriv_mg_tune;hmc_tm;offline_measurement")

include_directories(
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/lib/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<$<BOOL:${TM_USE_LEMON}>:${TM_LEMON_INCLUDE_DIRS}>
  ${TM_CLIME_INCLUDE_DIRS})

foreach(_prog ${tmlqcd_prog})
  add_executable(${_prog} "${_prog}.c")

  target_link_libraries(${_prog} PUBLIC hmc)
  set_target_properties(
    ${_prog}
    PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
               POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE "CXX")
endforeach()

if(TM_ENABLE_TESTS)
  list(
    APPEND
    tmlqcd_test_prog
    "check_locallity.c;hopping_test.cscalar_prod_r_test.c;test_eigenvalues.c;test_lemon.c"
  )
  if(TM_USE_LEMON)
    list(APPEND tmlqcd_test_prog test_lemon.c)
  endif()
  if(TN_USE_QPHIX)
    list(APPEND tmlqcd_test_prog qphix_test_Dslash.c)
  endif()

  foreach(_prog ${tmlqcd_test_prog})
    add_executable(${_prog} "${CMAKE_SOURCE_DIR}/src/bin/tests/${_prog}.c")

    target_link_libraries(${_prog} PUBLIC hmc)
    set_target_properties(
      ${_prog}
      PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                 POSITION_INDEPENDENT_CODE ON
                 LINKER_LANGUAGE "CXX")
  endforeach()
endif()
