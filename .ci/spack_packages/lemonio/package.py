# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems import cmake
from spack_repo.builtin.build_systems.cmake import CMakePackage, generator


from spack.package import *

class Lemonio(CMakePackage):
    """LEMON: Lightweight Parallel I/O library for Lattice QCD."""

    homepage = "https://github.com/etmc/lemon"
    git      = "https://github.com/etmc/lemon.git"
    license("GPL-3.0-or-later")

    version('master', branch='master')

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    depends_on('mpi')
    generator("ninja")

    def configure_args(self):
        args = []
        return args
