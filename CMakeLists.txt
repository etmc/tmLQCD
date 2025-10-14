cmake_minimum_required(VERSION 3.24)

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
include(CheckSymbolExists)
include(CheckLibraryExists)
include(CheckFunctionExists)

cmake_policy(SET CMP0048 NEW)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

project(
  tmlqcd
  DESCRIPTION "tmlQCD"
  HOMEPAGE_URL "http://www.itkp.uni-bonn.de/~urbach/software.html"
  VERSION "6.0.0"
  LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

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

# ##############################################################################
# Define the paths for static libraries and executables
# ##############################################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${cp2k_BINARY_DIR}/lib
  CACHE PATH "Single output directory for building all libraries.")

# Search for rocm in common locations
foreach(__var ROCM_ROOT CRAY_ROCM_ROOT ORNL_ROCM_ROOT CRAY_ROCM_PREFIX
    ROCM_PREFIX CRAY_ROCM_DIR)
  if($ENV{${__var}})
    list(APPEND CMAKE_PREFIX_PATH $ENV{__var})
    set(ROCM_PATH
      $ENV{__var}
      CACHE PATH "Path to ROCm installation")
  endif()
endforeach()

option(CMAKE_POSITION_INDEPENDENT_CODE "Enable position independent code" ON)

option(TMLQCD_USE_FFTW "Enable fftw support" OFF)
option(TMLQCD_USE_MPI "Enable MPI support" OFF)
option(TMLQCD_USE_CUDA "Enable QUDA support" OFF)
option(TMLQCD_USE_CUDA_HIP "Enable QUDA support" OFF)
option(TMLQCD_USE_HIP "Enable HIP support" OFF)
option(TMLQCD_USE_DDALPHAAMG "Enable DDalphaAMG support" OFF)
option(TMLQCD_USE_QPX "Enable QPX intrinsics on BGQ" OFF)
option(TMLQCD_USE_SPI "use IBM BG/Q SPI for communications" OFF)
option(TMLCDQ_USE_OPENMP "Enable openMP" ON)
option(TMLQCD_FIX_VOLUME "fix volume at compile time" OFF)
option(TMLQCD_ENABLE_ALIGNMENT "Automatically or expliclty align arrays to byte number. auto, none, 16, 32, 64" "auto")
option(TMLQCD_BGL_DRAM "use BGL dram window (BGL only!)" ON)
option(TMLQCD_USE_OPTIMIZATION "enable optimisation" ON)
option(TMLQCD_USE_GAUGE_COPY "Enable use of a copy of the gauge field" ON)
option(TMLQCD_USE_HALF_SPINOR "Use a Dirac Op. with halfspinor exchange" ON)
option(TMLQCD_USE_TSPLITPAR "Enable timeslice-splitted communications" ON)
option(TMLQCD_USE_LAPH "Enable computation of LapH eigensystem" OFF)
option(TMLQCD_USE_QPHIX "enable QPhiX" OFF)
option(TMLQCD_USE_SHMEM "Use shmem API" OFF)
option(TMLQCD_USE_INDEX_INDEP_GEOMETRY "enable Index independent addressing" OFF)
# asm instruction set
option(TMLQCD_USE_SSE2 "Use SSE2 instructions" OFF)
option(TMLQCD_USE_SSE3 "Use SSE3 instructions" OFF)
option(TMLQCD_USE_OPTERON "Use Opteron instructions" OFF)
option(TMLQCD_USE_SHMEM "Use shmem API" OFF)

# MPI dependent options
cmake_dependent_option(TMLQCD_PERSISTENT_MPI "Use persistent MPI calls for halfspinor [default=no]" OFF "TMLQCD_USE_MPI" OFF)
cmake_dependent_option(TMLQCD_NONBLOCKING_MPI "Use non-blocking MPI calls for spinor and gaug" ON "TMLQCD_USE_MPI" OFF)

# need to do it properly. Just a place holder
cmake_dependent_option(TMLQCD_MPI_DIMENSION "Use n dimensional parallelisation [default=4]" 4 "TMLQCD_USE_MPI" OFF)

# HIP dependent options
option_dependent_option(TMLQCD_USE_CUDA_HIP "Enable CUDA support in HIP" OFF "TMLQCD_USE_HIP" OFF)

# search for blas and lapack
find_package(BLAS REQUIRED)
#
find_package(LAPACK REQUIRED)

find_package(FLEX REQUIRED)
# do we need bison ?
find_package(BISON REQUIRED)

cmake_dependent_option(TMLQCD_USE_QUDA "Enable QUDA support" OFF "TMLQCD_USE_CUDA || TMLQCD_USE_HIP" OFF)
cmake_dependent_option(TMLQCD_USE_QUDA_EXPERIMENTAL "Enable QUDA support" OFF "TMLQCD_USE_QUDA" OFF)
cmake_dependent_option(TMLQCD_QUDA_FERMIONIC_FORCES "Enable support for fermionic forces using QUDA" ON "TMLQCD_USE_QUDA" OFF)

if (TMLQCD_USE_CUDA AND TMLQCD_USE_HIP)
  message(ERROR "HIP and CUDA are mutually exclusive. Please choose one GPU support only")
endif()

if (TMLQCD_USE_CUDA)
  enable_language(cuda)
endif()

if (TMLQCD_USE_HIP)
  enable_language(hip)

  # we may want to use hip-cuda for development or debugging purposes especially if AMD GPU is not possible. So allow it
  if (TMLQCD_USE_CUDA_HIP)
    find_package(CUDA)
  endif()

  if(CMAKE_HIP_PLATFORM MATCHES "amd")
  set(TMLQCD_GPU_PLATFORM_DFLAGS "__HIP_PLATFORM_AMD__")
  else()
  set(TMLQCD_GPU_PLATFORM_DFLAGS "__HIP_PLATFORM_NVIDIA__")
  endif()
endif()

if (TMLQCD_USE_QUDA)
  find_package(QUDA REQUIRED)
endif()

if (TMLQCD_USE_MPI)
  find_package(MPI REQUIRED)
endif()

if (TMLQCD_USE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C CXX)
endif()

if(TMLQCD_USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
endif()

# check for fftw3 (rely on pkgconfig).
if (TMLQCD_USE_FFTW)
  pkg_search_module(tmlqcd_fftw3 IMPORTED_TARGET GLOBAL fftw3)
  if (tmlqcd_fftw3_FOUND)
  add_library(tmlqcd::fftw3 ALIAS PkgConfig::tmlqcd_fftw3)
  endif()
endif()

# check for the presence of clock_gettime in libc or librt

check_symbol_exists(clock_gettime "time.h" HAVE_CLOCK_GETTIME)
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME_IN_RT)
check_function_exists(fseeko TMLQCD_HAVE_FSEEKO)
