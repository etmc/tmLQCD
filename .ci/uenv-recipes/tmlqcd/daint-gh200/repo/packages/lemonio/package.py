# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems import cmake
from spack_repo.builtin.build_systems.cmake import CMakePackage, generator


from spack.package import *

class Lemonio(AutotoolsPackage, CMakePackage):
    """LEMON: Lightweight Parallel I/O library for Lattice QCD."""

    homepage = "https://github.com/etmc/lemon"
    git      = "https://github.com/etmc/lemon.git"
    license("GPL-3.0-or-later")

    version('master', branch='master')

    depends_on("libtool", type="build", when="@master build_system=cmake")
    depends_on("cmake@4", type="build", when="master build_system=cmake")

    depends_on('mpi')

    generator("ninja")

class CMakeBuilder(cmake.CMakeBuilder):
    def cmake_args(self):
        spec = self.spec
        args = [
            self.define_from_variant("DBUILD_SHARED_LIBS", "shared"),
        ]
        return args

