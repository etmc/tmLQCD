# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems import cmake
from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack_repo.builtin.build_systems.rocm import ROCmPackage
from spack_repo.builtin.build_systems.cuda import CudaPackage

from spack.package import *


class Tmlqcd(CMakePackage, CudaPackage, ROCmPackage):
    """Base class for building tmlQCD."""

    homepage = "https://www.itkp.uni-bonn.de/~urbach/software.html"
    url = "https://github.com/etmc/tmLQCD/archive/refs/tags/rel-5-1-6.tar.gz"

    # todo: change this back to etmc as soon as cmake PR is merged
    git = "https://github.com/mtaillefumier/tmLQCD.git"
    license("GPL-3.0-or-later")

    maintainers("mtaillefumier")
    version("master", branch="master")

    # todo: remove this version as soon as
    # https://github.com/etmc/tmLQCD/pull/664 is merged
    version("cmake_support", branch="cmake_support")

    variant("lemon", default=False, description="Enable the lemon backend")
    variant("mpi", default=True, description="Enable mpi support")
    variant("DDalphaAMG", default=False, description="Enable DAlphaAMG support")
    variant("openmp", default=True, description="Enable OpenMP")
    variant("fftw", default=True, description="Enable FFTW interface")
    variant(
        "persistent_mpi",
        default=True,
        description="Enable persistent mpi calls for spinor and gauge fields",
        when="+mpi",
    )
    variant(
        "nonblocking_mpi",
        default=True,
        description="Enable non-blocking mpi calls for spinor and gauge fields",
        when="+mpi",
    )
    variant("fixedvolume", default=True, description="Enable fixed volume at compile time")
    variant(
        "alignment",
        default="auto",
        values=("none", "auto", "16", "32", "64"),
        description="Automatically or expliclty align arrays",
    )
    variant("gauge_copy", default=True, description="Enable gauge field copy")
    variant("half_spinor", default=True, description="Use a Dirac operator with half-spinor")
    variant("shared", default=False, description="Enable shared library")
    variant("shmem", default=False, description="Use shmem API")
    variant("quda", default=True, description="Enable the QUDA library", when="+cuda")
    variant("quda", default=True, description="Enable the QUDA library", when="+rocm")
    variant(
        "QPhiX", default=False, description="Enable the QPhiX library for Intel Xeon and Xeon Phis"
    )
    variant(
        "mpi_dimensions",
        default="4",
        values=("1", "2", "3", "4", "x", "xy", "xyz"),
        description="number of dimensions the mpi processes are distributed. the default is parallelization over all four dimensions txyz",
        when="+mpi",
    )

    generator("ninja")

    # language dependencies
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    # conflicts
    conflicts("+cuda", when="cuda_arch=none")
    conflicts("+rocm", when="amdgpu_target=none")
    conflicts("+cuda +rocm", msg="CUDA and ROCm support are mutually exclusive")

    # hard dependencies
    depends_on("c-lime")
    depends_on("blas")
    depends_on("lapack")
    depends_on("pkgconfig", type="build")

    # dependencies
    depends_on("mpi", when="+mpi")
    depends_on("lemonio", when="+lemon")

    depends_on("llvm-openmp", when="+rocm+openmp")

    with when("+quda"):
        depends_on(
            "quda+shared+twisted_mass+twisted_clover+clover+ndeg_twisted_clover+ndeg_twisted_mass+wilson+qdp+multigrid"
        )

        depends_on("quda+mpi", when="+mpi")
        depends_on("quda+cuda", when="+cuda")
        depends_on("quda+rocm", when="+rocm")

    depends_on("fftw-api@3", when="+fftw")


class CMakeBuilder(cmake.CMakeBuilder):
    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("TM_USE_LEMON", "lemon"),
            self.define_from_variant("TM_USE_MPI", "mpi"),
            self.define_from_variant("TM_USE_QUDA", "quda"),
            self.define_from_variant("TM_USE_CUDA", "cuda"),
            self.define_from_variant("TM_USE_HIP", "rocm"),
            self.define_from_variant("TM_USE_FFTW", "fftw"),
            self.define_from_variant("TM_USE_OMP", "openmp"),
            self.define_from_variant("TM_USE_SHMEM", "shmem"),
            self.define_from_variant("TM_USE_GAUGE_COPY", "gauge_copy"),
            self.define_from_variant("TM_USE_HALFSPINOR", "half_spinor"),
        ]

        # Use hipcc is case of a ROCm build
        if "+rocm" in self.spec:
            hip = self.spec["hip"]
            args.append(self.define("CMAKE_C_COMPILER", hip.hipcc))
            args.append(self.define("CMAKE_CXX_COMPILER", hip.hipcc))

            # help hipcc find openmp
            if "+openmp" in self.spec:
                omp = self.spec["llvm-openmp"]
                args.append(self.define("OpenMP_C_FLAGS", "-fopenmp"))
                args.append(self.define("OpenMP_CXX_FLAGS", "-fopenmp"))
                args.append(self.define("OpenMP_C_LIB_NAMES", "omp"))
                args.append(self.define("OpenMP_CXX_LIB_NAMES", "omp"))
                args.append(self.define("OpenMP_omp_LIBRARY", "{0}/libomp.so".format(omp.prefix.lib)))
                args.append(self.define("OpenMP_CXX_INCLUDE_DIR", omp.prefix.include))

        return args
