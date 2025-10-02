# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os

from spack_repo.builtin.build_systems.autotools import AutotoolsPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage
from spack_repo.builtin.build_systems.cuda import CudaPackage

import llnl.util.lang

from spack.package import *
class Tmlqcd(AutotoolsPackage, CudaPackage, ROCmPackage):
    """Base class for building tmlQCD."""

    homepage = "https://www.itkp.uni-bonn.de/~urbach/software.html"
    url = "https://github.com/etmc/tmLQCD/archive/refs/tags/rel-5-1-6.tar.gz"
    git = "https://github.com/etmc/tmLQCD.git"
    license("GPL-3.0-or-later")

    maintainers("mtaillefumier")

    version("5-1-6", sha256="9af2930cae47acb2f8b5155525590ecda4cc8ececac42845e94617c75a9a9e19")
    version("master", branch="master")

    variant("lime", default=True, description="Enable the lime backend")
    variant("mpi", default=True, description="Enable mpi support")
    variant("DDalphaAMG", default=False, description="Enable DAlphaAMG support")
    variant("qpx", default=False, description="Enable IBM qpx intrinsic")
    variant("spi", default=False, description="Enable IBM SPI communication")
    variant("openmp", default=True, description="Enable OpenMP")
    variant("fftw", default=True, description="Enable FFTW interface")
    variant(
        "persistentmpi",
        default=True,
        description="Enable persistent mpi calls for spinor and gauge fields",
        when="+mpi",
    )
    variant(
        "nonblockingmpi",
        default=True,
        description="Enable non-blocking mpi calls for spinor and gauge fields",
        when="+mpi",
    )
    variant("fixedvolume", default=True, description="Enable fixed volume at compile time")
    variant("kojak", default=False, description="instrumentalise for KOJAK")
    variant(
        "alignment",
        default="auto",
        values=("none", "auto", "16", "32", "64"),
        description="Automatically or expliclty align arrays",
    )
    variant("gauge_copy", default=True, description="Enable gauge field copy")
    variant("half_spinor", default=True, description="Use a Dirac operator with half-spinor")
    variant("shmem", default=False, description="Use shmem API")
    variant("indexindepgeo", default=False, description="enable Index independent addressing")
    variant("tsplitpar", default=False, description="Enable timeslice-splitted communications")
    variant("laph", default=False, description="enable computation of LapH eigensystem")
    variant("quda", default=True, description="Enable the QUDA library", when="+cuda",)
    variant("quda", default=True, description="Enable the QUDA library", when="+rocm",)
    variant(
        "quda_experimental",
        default=True,
        description="Enable QUDA experimental features",
        when="+quda",
    )
    variant(
        "quda_fermionic_forces",
        default=True,
        description="Enable support for fermionic forces using QUDA",
        when="+quda",
    )
    variant(
        "QPhiX", default=False, description="Enable the QPhiX library for Intel Xeon and Xeon Phis"
    )
    variant(
        "mpi_dimensions",
        default="4",
        values=("1", "2", "3", "4", "x", "xy", "xyz"),
        description="number of dimentions the mpi processes are distributed",
        when="+mpi",
    )
    
    variant("lemon", default=False, description="Enable the lemon backend", when="+mpi")
    # language dependencies
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    # conflicts
    conflicts("+cuda", when="cuda_arch=none")
    conflicts("+rocm", when="amdgpu_target=none")

    # dependencies
    depends_on("c-lime", when="lime")
    depends_on("mpi", when="+mpi")

    with when("+quda"):
        depends_on(
            "quda+twisted_mass+twisted_clover+clover+ndeg_twisted_clover+ndeg_twisted_mass+wilson+qdp+staggered+usqcd+qio+qmp"
        )

        depends_on("quda+mpi", when="+mpi")
        depends_on("quda+cuda", when="+cuda")
        depends_on("quda+rocm", when="+rocm")

    depends_on("fftw-api@3", when="+fftw")
    depends_on("blas")
    depends_on("lapack")
    depends_on("pkgconfig", type="build")
    depends_on("etmc-lemon", when="lemon")

    @property
    def force_autoreconf(self):
        return

    def configure_args(self):
        configure_args = []
        configure_args.append("--with-lapack=%s" % self.spec["blas"].libs)

        if self.spec.satisfies("+lemon"):
            configure_args.append(" --with-lemon=%s" % self.spec["lemon"].prefix)

        if self.spec.satisfies("+cuda"):
            configure_args.append(" --with-cudadir=%s" % self.spec["cuda"].prefix)
            configure_args.extend(self.enable_or_disable("cuda_experimental"))

        if self.satisfies("+rocm"):
            configure_args.append(" --with-hipdir=%s" % self.spec["rocm"].prefix)

        configure_args.extend(self.enable_or_disable("mpi"))
        configure_args.extend(self.enable_or_disable("sse2"))
        configure_args.extend(self.enable_or_disable("sse3"))
        configure_args.extend(self.enable_or_disable("openmp"))
        configure_args.extend(self.enable_or_disable("persistentmpi"))
        configure_args.extend(self.enable_or_disable("nonblockingmpi"))
        configure_args.extend(self.enable_or_disable("quda_fermionic_forces"))
        configure_args.extend(self.enable_or_disable("gaugecopy"))
        configure_args.extend(self.enable_or_disable("qpx"))
        configure_args.extend(self.enable_or_disable("spi"))
        configure_args.extend(self.enable_or_disable("fftw"))
        configure_args.extend(self.enable_or_disable("optimize"))
        configure_args.extend(self.enable_or_disable("shared"))
        configure_args.extend(self.with_or_without("pic"))
        configure_args.extend(self.enable_or_disable("indexindepgeo"))
        configure_args.extend(self.enable_or_disable("halfspinor"))
        configure_args.extend(self.enable_or_disable("tsplitpar"))
        configure_args.extend(self.with_or_without("fixedvolume"))
        if self.spec.satisfies("+mpi"):
            configure_args.append("--with-mpidimension=%s" % self.spec["mpi_dimensions"])

        configure_args.append("--enable-alignment=%s" % self.spec["alignment"].value)
        return configure_args
