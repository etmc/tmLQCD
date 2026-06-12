list(APPEND tmlqcd_prog
     "invert;benchmark;deriv_mg_tune;hmc_tm;offline_measurement")

include_directories(
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/lib/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<$<BOOL:${TM_USE_LEMON}>:${TM_LEMON_INCLUDE_DIRS}>
  ${TM_CLIME_INCLUDE_DIRS})

foreach(_prog ${tmlqcd_prog})
  add_executable(${_prog} "${_prog}.c")

  target_link_libraries(${_prog} PUBLIC tmlqcd)
  set_target_properties(
    ${_prog}
    PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
               POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE "CXX")
endforeach()

# check_locallity; compilation failure scalar_prod_r_test: compilation error
# header file missing test_eigenvalues: observables.h missing
if(TM_ENABLE_TESTS)
  list(APPEND tmlqcd_test_prog "hopping_test")
  if(TM_USE_LEMON)
    list(APPEND tmlqcd_test_prog "test_lemon")
  endif()
  if(TN_USE_QPHIX)
    list(APPEND tmlqcd_test_prog "qphix_test_Dslash")
  endif()

  foreach(_prog ${tmlqcd_test_prog})
    add_executable(${_prog} "${CMAKE_SOURCE_DIR}/src/bin/tests/${_prog}.c")
    target_include_directories(${_prog}
                               PUBLIC "${CMAKE_SOURCE_DIR}/src/lib/test")
    target_link_libraries(${_prog} PUBLIC tmlqcd)
    set_target_properties(
      ${_prog}
      PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                 POSITION_INDEPENDENT_CODE ON
                 LINKER_LANGUAGE "CXX")
  endforeach()
endif()

install(TARGETS ${tmlqcd_prog} RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
