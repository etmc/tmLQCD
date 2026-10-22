list(APPEND tmlqcd_tests
     "test_buffers;test_clover;test_linalg;test_rat;test_sample;test_su3")

if(TM_ENABLE_TESTS)
  add_executable(
    test_buffers
    "test_buffers.c;test_buffers_gauge.c;${CMAKE_SOURCE_DIR}/src/lib/cu/cu.c")
  add_executable(
    test_clover
    "test_clover.c;test_clover_six_invert.c;${CMAKE_SOURCE_DIR}/src/lib/cu/cu.c"
  )
  add_executable(
    test_linalg
    "test_linalg.c;test_linalg_spinor.c;${CMAKE_SOURCE_DIR}/src/lib/cu/cu.c")
  add_executable(
    test_rat "test_rat.c;test_rat_init.c;${CMAKE_SOURCE_DIR}/src/lib/cu/cu.c")
  add_executable(
    test_sample
    "test_sample.c;test_sample_ts1.c;test_sample_ts2.c;${CMAKE_SOURCE_DIR}/src/lib/cu/cu.c"
  )
  add_executable(
    test_su3
    "test_su3_algebra.c;test_su3.c;${CMAKE_SOURCE_DIR}/src/lib/cu/cu.c")

  foreach(_prog ${tmlqcd_tests})
    target_link_libraries(${_prog} PUBLIC tmlqcd)
    target_include_directories(
      ${_prog}
      PUBLIC "${CMAKE_SOURCE_DIR}/src/lib/test" ${CMAKE_BINARY_DIR}
             ${CMAKE_SOURCE_DIR}/src/lib ${CMAKE_SOURCE_DIR}/src/lib/include
             ${CMAKE_SOURCE_DIR}/src/lib/cu)
    set_target_properties(
      ${_prog}
      PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                 POSITION_INDEPENDENT_CODE ON
                 LINKER_LANGUAGE "CXX")
  endforeach()
endif()
