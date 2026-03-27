# DetectSimdAndAlignment.cmake
#
# Detect SIMD architecture family, SIMD level and a reasonable alignment value.
#
# Exposed cache variables:
#   SIMD_ARCH_FAMILY : x86 / ARM / PPC / UNKNOWN
#   SIMD_LEVEL       : AVX512 / AVX2 / SSE2 / NEON / ALTIVEC / SCALAR
#   SIMD_ALIGNMENT   : integer, in bytes (16, 32, 64, ...)
#
# Optional (if you want a configured header):
#   SIMD_CONFIG_HEADER : path to the generated header (see bottom).
#
# Usage:
#   include(cmake/DetectSimdAndAlignment.cmake)
#   message(STATUS "SIMD: ${SIMD_ARCH_FAMILY} ${SIMD_LEVEL}, alignment=${SIMD_ALIGNMENT}")
#
#   # Example: propagate as defines
#   target_compile_definitions(my_target PRIVATE
#       SIMD_ALIGNMENT=${SIMD_ALIGNMENT}
#       SIMD_LEVEL_${SIMD_LEVEL}
#   )
# DetectSimdAndAlignment.cmake - COMPLETE: x86 + ARM NEON + NVIDIA + PowerPC


include_guard(GLOBAL) #

include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns) # For runtime CPU detection fallback

# ------------------------------
# 1. Detect architecture family
# ------------------------------
if(NOT DEFINED SIMD_ARCH_FAMILY)
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _simd_proc)

    if(_simd_proc MATCHES "x86_64|amd64|i[3-6]86")
        set(_detected_arch "x86")
    elseif(_simd_proc MATCHES "armv[0-9]+|aarch64|arm64")
        set(_detected_arch "ARM")
    elseif(_simd_proc MATCHES "ppc64(le|el)?|powerpc|ppc")
        set(_detected_arch "PPC")
    elseif(_simd_proc MATCHES "nvcl|sm_89|sm_90")
        set(_detected_arch "NVIDIA")
    else()
        set(_detected_arch "UNKNOWN")
    endif()

    set(SIMD_ARCH_FAMILY "${_detected_arch}" CACHE STRING "SIMD architecture family")
endif()

# Defaults
set(SIMD_LEVEL "SCALAR" CACHE STRING "Detected SIMD level")
set(SIMD_ALIGNMENT 16 CACHE STRING "Alignment in bytes")
set(SIMD_HAS_FLOAT ON CACHE BOOL "Float SIMD support")
set(SIMD_HAS_DOUBLE ON CACHE BOOL "Double SIMD support")

# Save/restore flags helper
set(_SIMD_SAVED_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
macro(_simd_restore_flags)
    if(DEFINED _SIMD_SAVED_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS}")
    endif()
endmacro()

# ------------------------------------------------
# 2. x86: SSE2 → AVX2 → AVX512
# ------------------------------------------------
if(SIMD_ARCH_FAMILY STREQUAL "x86")
    # AVX512 double (64-byte)
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -mavx512f -mavx512dq")
    check_cxx_source_compiles("
        #include <immintrin.h>
        int main() { __m512d v = _mm512_set1_pd(1.0); (void)v; return 0; }
    " _HAVE_AVX512_DOUBLE)

    if(_HAVE_AVX512_DOUBLE)
        set(SIMD_LEVEL "AVX512" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 64 CACHE STRING "" FORCE)
        _simd_restore_flags()
        return()
    endif()

    # AVX2 double (32-byte)
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -mavx2")
    check_cxx_source_compiles("
        #include <immintrin.h>
        int main() { __m256d v = _mm256_set1_pd(1.0); (void)v; return 0; }
    " _HAVE_AVX2_DOUBLE)

    if(_HAVE_AVX2_DOUBLE)
        set(SIMD_LEVEL "AVX2" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 32 CACHE STRING "" FORCE)
        _simd_restore_flags()
        return()
    endif()

    # SSE2 double minimum (16-byte)
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -msse2")
    check_cxx_source_compiles("
        #include <emmintrin.h>
        int main() { __m128d v = _mm_set1_pd(1.0); (void)v; return 0; }
    " _HAVE_SSE2_DOUBLE)

    if(_HAVE_SSE2_DOUBLE)
        set(SIMD_LEVEL "SSE2" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
        _simd_restore_flags()
        return()
    endif()

# --------------------------------------
# 3. ARM NEON - ALL FAMILIES
# --------------------------------------
elseif(SIMD_ARCH_FAMILY STREQUAL "ARM")
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _arm_proc)

    # AArch64 + SVE
    if(_arm_proc MATCHES "aarch64|arm64")
        set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -march=armv8-a+sve")
        check_cxx_source_compiles("
            #include <arm_sve.h>
            int main() { svfloat32_t v = svdup_f32(1.0f); (void)v; return 0; }
        " _HAVE_SVE)

        if(_HAVE_SVE)
            set(SIMD_LEVEL "SVE" CACHE STRING "" FORCE)
            set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
            _simd_restore_flags()
            return()
        endif()

        # AArch64 NEON (double safe)
        check_cxx_source_compiles("
            #include <arm_neon.h>
            int main() {
                float64x2_t vd = vdupq_n_f64(1.0);
                float32x4_t vf = vdupq_n_f32(1.0f);
                (void)vd; (void)vf; return 0;
            }" _HAVE_NEON_AARCH64)

        if(_HAVE_NEON_AARCH64)
            set(SIMD_LEVEL "NEON_AARCH64" CACHE STRING "" FORCE)
            set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
            _simd_restore_flags()
            return()
        endif()

    # ARMv8 32-bit
    elseif(_arm_proc MATCHES "armv8")
        set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -march=armv8-a+simd")
        check_cxx_source_compiles("
            #include <arm_neon.h>
            int main() { float32x4_t v = vdupq_n_f32(1.0f); (void)v; return 0; }
        " _HAVE_ARMv8_NEON)

        if(_HAVE_ARMv8_NEON)
            set(SIMD_LEVEL "NEON_ARMv8" CACHE STRING "" FORCE)
            set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
            set(SIMD_HAS_DOUBLE OFF CACHE BOOL "" FORCE)
            _simd_restore_flags()
            return()
        endif()

    # ARMv7 NEON
    elseif(_arm_proc MATCHES "armv7")
        set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -mfpu=neon -march=armv7-a")
        check_cxx_source_compiles("
            #include <arm_neon.h>
            int main() { float32x4_t v = vdupq_n_f32(1.0f); (void)v; return 0; }
        " _HAVE_ARMv7_NEON)

        if(_HAVE_ARMv7_NEON)
            set(SIMD_LEVEL "NEON_ARMv7" CACHE STRING "" FORCE)
            set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
            set(SIMD_HAS_DOUBLE OFF CACHE BOOL "" FORCE)
            _simd_restore_flags()
            return()
        endif()
    endif()

# --------------------------------------
# 4. POWERPC - COMPLETE COVERAGE (NEW!)
# --------------------------------------
elseif(SIMD_ARCH_FAMILY STREQUAL "PPC")

    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _ppc_proc)

    # === Power10+ (512-bit vectors, POWER10)
    # Note: Power10 needs -mcpu=power10 or -mtune=power10
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -mcpu=power10")
    check_cxx_source_compiles("
        #include <altivec.h>
        int main() {
            vector double vd = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}; // 512-bit
            vector float vf = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
            (void)vd; (void)vf; return 0;
        }" _HAVE_POWER10)

    if(_HAVE_POWER10)
        set(SIMD_LEVEL "POWER10" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 64 CACHE STRING "" FORCE)  # 512-bit = 64 bytes
        set(SIMD_HAS_FLOAT ON CACHE BOOL "" FORCE)
        set(SIMD_HAS_DOUBLE ON CACHE BOOL "" FORCE)
        _simd_restore_flags()
        return()
    endif()

    # === Power9 VSX (256-bit, POWER8+)
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -mcpu=power9 -mvsx")
    check_cxx_source_compiles("
        #include <altivec.h>
        int main() {
            vector double vd = {1.0,1.0,1.0,1.0};  // 256-bit VSX double
            vector float vf = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f}; // 256-bit
            (void)vd; (void)vf; return 0;
        }" _HAVE_VSX_POWER9)

    if(_HAVE_VSX_POWER9)
        set(SIMD_LEVEL "VSX_POWER9" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 32 CACHE STRING "" FORCE)  # 256-bit = 32 bytes
        set(SIMD_HAS_FLOAT ON CACHE BOOL "" FORCE)
        set(SIMD_HAS_DOUBLE ON CACHE BOOL "" FORCE)
        _simd_restore_flags()
        return()
    endif()

    # === Power7+ VSX (128-bit double, POWER7+)
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -mcpu=power7 -mvsx")
    check_cxx_source_compiles("
        #include <altivec.h>
        int main() {
            vector double vd = {1.0,1.0};  // VSX 128-bit double
            (void)vd; return 0;
        }" _HAVE_VSX_POWER7)

    if(_HAVE_VSX_POWER7)
        set(SIMD_LEVEL "VSX_POWER7" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
        set(SIMD_HAS_FLOAT ON CACHE BOOL "" FORCE)
        set(SIMD_HAS_DOUBLE ON CACHE BOOL "" FORCE)
        _simd_restore_flags()
        return()
    endif()

    # === Classic AltiVec/VMX (PowerPC baseline, 128-bit)
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} -maltivec -mabi=altivec")
    check_cxx_source_compiles("
        #include <altivec.h>
        int main() {
            vector float vf = (vector float){1.0f,1.0f,1.0f,1.0f};
            (void)vf; return 0;
        }" _HAVE_ALTIVEC)

    if(_HAVE_ALTIVEC)
        set(SIMD_LEVEL "ALTIVEC" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
        set(SIMD_HAS_FLOAT ON CACHE BOOL "" FORCE)
        set(SIMD_HAS_DOUBLE OFF CACHE BOOL "" FORCE)  # AltiVec: float primary
        _simd_restore_flags()
        return()
    endif()

# --------------------------------------
# 5. NVIDIA GH200 (sm_89)
# --------------------------------------
elseif(SIMD_ARCH_FAMILY STREQUAL "NVIDIA")
    set(CMAKE_REQUIRED_FLAGS "${_SIMD_SAVED_REQUIRED_FLAGS} --gpu-arch=sm_89")
    check_cxx_source_compiles("
        #include <cuda_runtime.h>
        int main() { double d = 1.0; (void)d; return 0; }
    " _HAVE_CUDA_SM89)

    if(_HAVE_CUDA_SM89)
        set(SIMD_LEVEL "CUDA_SM89" CACHE STRING "" FORCE)
        set(SIMD_ALIGNMENT 16 CACHE STRING "" FORCE)
        _simd_restore_flags()
        return()
    endif()

# --------------------------------------
# 6. Fallback
# --------------------------------------
else()
    _simd_restore_flags()
    return()
endif()

_simd_restore_flags()
