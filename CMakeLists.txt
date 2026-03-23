cmake_minimum_required(VERSION 3.30)

project(
  tmlqcd
  DESCRIPTION "tmlQCD"
  HOMEPAGE_URL "https://github.com/etmc/tmLQCD"
  VERSION "6.0.0"
  LANGUAGES C CXX)

# include our cmake snippets
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# =================================================================================================
# REQUIRE OUT-OF-SOURCE BUILDS
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(
    FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory."
  )
endif()

# =================================================================================================
# PROJECT AND VERSION
include(CMakeDependentOption)
include(GNUInstallDirs)
include(FetchContent)
include(CMakePackageConfigHelpers)

cmake_policy(SET CMP0048 NEW)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
  set(CMAKE_C_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_HIP_STANDARD)
  set(CMAKE_HIP_STANDARD 14)
  set(CMAKE_HIP_STANDARD_REQUIRED ON)
endif()

find_package(PkgConfig)

# Search for rocm in common locations
foreach(
  __var
  ROCM_ROOT
  HIP_ROOT
  HIP_PATH
  CRAY_ROCM_ROOT
  ORNL_ROCM_ROOT
  CRAY_ROCM_PREFIX
  ROCM_PREFIX
  CRAY_ROCM_DIR
  ROCM_PATH)
  if($ENV{${__var}})
    list(APPEND CMAKE_PREFIX_PATH $ENV{__var})
    set(ROCM_PATH
        $ENV{__var}
        CACHE PATH "Path to ROCm installation")
  endif()
endforeach()

option(CMAKE_POSITION_INDEPENDENT_CODE "Enable position independent code" ON)
option(BUILD_SHARED_LIBS "Enable shared library" OFF)
option(TM_USE_FFTW "Enable fftw support" OFF)
option(TM_USE_MPI "Enable MPI support" OFF)
option(TM_USE_CUDA "Enable QUDA support" OFF)
option(TM_USE_HIP "Enable HIP support" OFF)
option(TM_USE_DDalphaAMG "Enable DDalphaAMG support" OFF)
option(TM_USE_OMP "Enable OpenMP" ON)
option(TM_FIXEDVOLUME "Fix volume at compile time" OFF)
set(TM_ENABLE_ALIGNMENT
    "auto"
    CACHE
      STRING
      "Automatically or expliclty align arrays to byte number. auto, none, 16, 32, 64"
)

set_property(CACHE TM_ENABLE_ALIGNMENT PROPERTY STRINGS "auto" "none" "16" "32"
                                                "64")

option(TM_USE_OPTIMIZATION "enable optimisation" ON)
option(TM_USE_GAUGE_COPY "Enable use of a copy of the gauge field" ON)
option(TM_USE_HALFSPINOR "Use a Dirac Op. with halfspinor exchange" ON)
option(TM_USE_QPHIX "Enable QPhiX" OFF)
option(TM_USE_SHMEM "Use shmem API" OFF)
option(TM_USE_QUDA "Enable QUDA support" OFF)
option(TM_ENABLE_WARNINGS "Enable all warnings" ON)
option(TM_ENABLE_TESTS "Enable tests" OFF)
set(TM_QPHIX_SOALEN
    "4"
    CACHE STRING "QPhiX specific parameter")
# MPI dependent options
cmake_dependent_option(
  TM_PERSISTENT_MPI "Use persistent MPI calls for halfspinor [default=no]" OFF
  "TM_USE_MPI" OFF)
cmake_dependent_option(
  TM_NONBLOCKING_MPI "Use non-blocking MPI calls for spinor and gauge" ON
  "TM_USE_MPI" OFF)

# need to do it properly. Just a place holder
cmake_dependent_option(
  TM_MPI_DIMENSION "Use n dimensional parallelisation [default=4]" 4
  "TM_USE_MPI" OFF)

# HIP dependent options
cmake_dependent_option(TM_USE_CUDA_HIP "Enable CUDA support in HIP" OFF
                       "TM_USE_HIP" OFF)

# clime and lemon depend on MPI
cmake_dependent_option(TM_USE_LEMON "Use the lemon io library" OFF "TM_USE_MPI"
                       ON)

# GPU dependent options
cmake_dependent_option(TM_USE_NVHPC "Enable Nvidia HPC toolkit" OFF
                       "TM_USE_CUDA" OFF)

# DDAlphaAMG specific options

cmake_dependent_option(
  DDalphaAMG_ENABLE_PARAMOUNT_OUTPUT "Enable paramount output support" ON
  "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(DDalphaAMG_ENABLE_FGMRES_RESTEST "Enable GMRES test" OFF
                       "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(
  DDalphaAMG_ENABLE_PROFILING "Enable paramount output support" OFF
  "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(DDalphaAMG_ENABLE_TRACK_RES "Enable track res support"
                       ON "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(
  DDalphaAMG_ENABLE_SINGLE_ALLREDUCE_ARNOLDI OFF
  "Enable paramount output support" OFF "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(
  DDalphaAMG_ENABLE_COARSE_RES "Enable paramount output support" OFF
  "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(
  DDalphaAMG_ENABLE_SCHWARZ_RES "Enable paramount output support" OFF
  "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(DDalphaAMG_ENABLE_OMP "Enable OpenMP support" ON
                       "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(
  DDalphaAMG_ENABLE_TESTVECTOR_ANALYSIS "Enable vector analysis support" OFF
  "TM_USE_DDalphaAMG" OFF)

cmake_dependent_option(DDalphaAMG_ENABLE_HDF5 "Enable HDF5 support" OFF
                       "TM_USE_DDalphaAMG" OFF)

# search for blas and lapack
find_package(BLAS REQUIRED)
#
find_package(LAPACK REQUIRED)
set(TM_LAPACK ON)
find_package(FLEX REQUIRED)

set(PACKAGE_NAME ${PROJECT_DESCRIPTION})
set(PACKAGE_VERSION ${PROJECT_VERSION})
set(PACKAGE_TARNAME "tmlqcd")
set(PACKAGE_BUGREPORT "curbach@gmx.de")
set(PACKAGE_STRING "${PROJECT_DESCRIPTION} ${PROJECT_VERSION}")

set(ALIGN " ")
set(ALIGN_BASE "0")
set(ALIGN_BASE32 "0")
set(ALIGN32 " ")

# DO NOT MERGE the two if statements as otherwise the automatic alignment will
# not be taken into account

if(${TM_ENABLE_ALIGNMENT} STREQUAL "auto")
  include(cmake/DetectSimdAndAlignment.cmake)
  message(
    STATUS "SIMD: ${SIMD_LEVEL} (${SIMD_ARCH_FAMILY}), align=${SIMD_ALIGNMENT}")
endif()

if(${TM_ENABLE_ALIGNMENT} STREQUAL "none")
  set(ALIGN_BASE "0x00")
  set(ALIGN " ")
  set(ALIGN_BASE32 "0x00")
  set(ALIGN32 " ")
elseif((${TM_ENABLE_ALIGNMENT} STREQUAL "16") OR (${SIMD_ALIGNMENT} EQUAL 16))
  set(ALIGN_BASE "0x0F")
  set(ALIGN "__attribute__ ((aligned (16)))")
  set(ALIGN_BASE32 "0x0F")
  set(ALIGN32 "__attribute__ ((aligned (16)))")
elseif((${TM_ENABLE_ALIGNMENT} STREQUAL "32") OR (${SIMD_ALIGNMENT} EQUAL 32))
  set(ALIGN_BASE "0x2F")
  set(ALIGN "__attribute__ ((aligned (32)))")
  set(ALIGN_BASE32 "0x2F")
  set(ALIGN32 "__attribute__ ((aligned (32)))")
elseif((${TM_ENABLE_ALIGNMENT} STREQUAL "64") OR (${SIMD_ALIGNMENT} EQUAL 64))
  set(ALIGN_BASE "0x3F")
  set(ALIGN "__attribute__ ((aligned (64)))")
  set(ALIGN_BASE32 "0x3F")
  set(ALIGN32 "__attribute__ ((aligned (64)))")
else()
  message(
    FATAL_ERROR
      "Unusable value for array alignment. Allowed values are: auto, none, 16, 32, 64"
  )
endif()

if(TM_USE_MPI)
  find_package(MPI REQUIRED)
  if(TM_NONBLOCKING_MPI)
    set(TM_NONBLOCKING ON)
  endif()
  if(TM_PERSISTENT_MPI)
    set(TM_PERSISTENT ON)
  endif()
endif()

if(TM_USE_OMP)
  find_package(OpenMP REQUIRED COMPONENTS C CXX)
endif()

if(TM_USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
endif()

if(TM_USE_LEMON)
  FetchContent_Declare(
    lemon
    GIT_REPOSITORY https://github.com/etmc/lemon
    GIT_TAG 42c2d99dae9c04ca1f09c532a8f9bcb1bb667528
    FIND_PACKAGE_ARGS NAMES lemon)
  FetchContent_MakeAvailable(lemon)
endif()

find_package(CLime REQUIRED)
set(TM_USE_LIME ON)

if(TM_USE_QUDA)
  find_package(QUDA REQUIRED config)
endif()

if(TM_USE_SHMEM)
  message(INFO "SHMEM needs to be included")
endif()

if(TM_USE_CUDA AND TM_USE_HIP)
  message(
    ERROR
    "HIP and CUDA are mutually exclusive. Please choose one GPU support only")
endif()

if(TM_USE_CUDA OR QUDA_TARGET_CUDA)
  enable_language(CUDA)

  # placeholder for nvhpc for future use
  if(TM_USE_NVHPC)
    find_package(NVHPC REQUIRED COMPONENTS CUDA MATH HOSTUTILS NCCL)
  else()
    find_package(CUDAToolkit REQUIRED)
  endif()
endif()

# We may want to use hip-cuda for development or debugging purposes especially
# if AMD GPU access is not possible. So allow it

if(TM_USE_HIP OR QUDA_TARGET_HIP)
  enable_language(HIP)
  find_package(hipblas)
  find_package(hipfft)
  if(TM_USE_CUDA_HIP)
    find_package(CUDA)
  endif()

  if(CMAKE_HIP_PLATFORM MATCHES "amd")
    set(TM_GPU_PLATFORM_DFLAGS "__HIP_PLATFORM_AMD__")
  else()
    set(TM_GPU_PLATFORM_DFLAGS "__HIP_PLATFORM_NVIDIA__")
  endif()
endif()

if(TM_USE_QPHIX)
  find_package(Qphix REQUIRED)
endif()

# check for fftw3 (rely on pkgconfig).
if(TM_USE_FFTW)
  pkg_search_module(tmlqcd_fftw3 IMPORTED_TARGET GLOBAL fftw3)
  if(tmlqcd_fftw3_FOUND)
    add_library(tmlqcd::fftw3 ALIAS PkgConfig::tmlqcd_fftw3)
  endif()
endif()

if(TM_USE_DDalphaAMG)
  add_subdirectory(DDalphaAMG)
endif()

if(TM_ENABLE_WARNINGS)
  add_compile_options($<$<COMPILE_LANG_AND_ID:C,GNU>:-Wall>
                      $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wall>)
endif()

# set the parallelization

if(TM_USE_MPI)
  if(TM_MPI_DIMENSION STREQUAL "1")
    # T parallelisation
    set(TM_PARALLELT ON)
  elseif(TM_MPI_DIMENSION STREQUAL "2")
    # XT parallelisation
    set(TM_PARALLELXT ON)
  elseif(TM_MPI_DIMENSION STREQUAL "3")
    set(TM_PARALLELXYT ON)
    # XYZ parallelisation
  elseif(TM_MPI_DIMENSION STREQUAL "4")
    # timeslice-splitted communications
    set(TM_PARALLELXYZT ON)
  elseif(TM_MPI_DIMENSION STREQUAL "X")
    set(TM_PARALLELX ON)
  elseif(TM_MPI_DIMENSION STREQUAL "XY")
    set(TM_PARALLELXY ON)
  elseif(TM_MPI_DIMENSION STREQUAL "XYZ")
    set(TM_PARALLELXYZ ON)
  else()
    set(TM_PARALLELXYZT ON)
  endif()
endif()

if(TM_USE_HALFSPINOR AND NOT TM_USE_GAUGE_COPY)
  message(
    FATAL_ERROR
      "The TM_USE_GAUGE_COPY option should also be set to ON when TM_USE_HALFSPINOR is ON"
  )
endif()
# keep the autotool config.h header.
configure_file("${PROJECT_SOURCE_DIR}/cmake/tmlqcd_config_internal.h.in"
               "${PROJECT_BINARY_DIR}/tmlqcd_config_internal.h" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/lib/fixed_volume.h.in"
               "${PROJECT_BINARY_DIR}/fixed_volume.h" @ONLY)
# check if git command exists
find_program(GIT_EXE NAMES git)

# generate version header
string(TIMESTAMP TM_TIMESTAMP "%Y-%m-%d %H:%M:%S")
if(DEFINED GIT_EXE AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse HEAD
    OUTPUT_VARIABLE TM_SHA
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git describe --all
    OUTPUT_VARIABLE TM_GIT_BRANCH
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "git hash ${TM_SHA}")
else()
  set(TM_SHA
      "https://github.com/etmc/tmLQCD/releases/tag/rel-${TMLQCD_VERSION_MAJOR}-${TMLQCD_VERSION_MINOR}"
  )
endif()

configure_file(cmake/git_hash.c.in git_hash.c @ONLY)
add_subdirectory(src/lib)
add_subdirectory(src/bin)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/tmlQCDConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)

configure_file("${PROJECT_SOURCE_DIR}/cmake/tmlQCD-config.cmake.in"
               "${PROJECT_BINARY_DIR}/tmlQCD-config.cmake" @ONLY)

install(FILES "${PROJECT_BINARY_DIR}/tmlQCD-config.cmake"
              "${PROJECT_BINARY_DIR}/tmlQCDConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/tmlQCD")

install(FILES "${PROJECT_BINARY_DIR}/libtmlQCD.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/tmlQCD"
  FILES_MATCHING
  PATTERN "*.cmake")
