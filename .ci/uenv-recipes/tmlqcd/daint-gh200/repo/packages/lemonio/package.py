# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.autotools import AutotoolsPackage


from spack.package import *

class Lemonio(AutotoolsPackage):
    """LEMON: Lightweight Parallel I/O library for Lattice QCD."""

    homepage = "https://github.com/etmc/lemon"
    git      = "https://github.com/etmc/lemon.git"
    license("GPL-3.0-or-later")

    version('master', branch='master')

    depends_on("autoconf", type="build", when="@master build_system=autotools")
    depends_on("automake", type="build", when="@master build_system=autotools")
    depends_on("libtool", type="build", when="@master build_system=autotools")

    depends_on('mpi')

    def configure_args(self):
        args = []
        args.append('CC={0}'.format(self.spec['mpi'].mpicc))
        return args
