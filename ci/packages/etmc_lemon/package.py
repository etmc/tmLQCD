# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.autotools import AutotoolsPackage

from spack.package import *

class EtmcLemon(AutotoolsPackage):
    """Base class for building tmlQCD."""

    homepage = "https://github.com/etmc/lemon"
    url = "https://github.com/etmc/lemon/archive/refs/tags/r1.1.tar.gz"
    git = "https://github.com/etmc/lemon.git"
    license("GPL-3.0-or-later")

    maintainers("mtaillefumier")

    version("1.1", sha256="c74cc458b0e17bed81796b1e6b277c8f4e13b81003ff9af9791ff41bc18b46b6")
    version("master", branch="master")
    
    # language dependencies
    depends_on("c", type="build")

    # build system dependencies
    depends_on("autoconf", type="build")
    depends_on("automake", type="build")
    depends_on("libtool", type="build")
    depends_on("m4", type="build")

    # dependencies
    depends_on("mpi")

    def configure_args(self):
        configure_args = [
                "CC={0}".format(self.spec["mpi"].mpicc)]

        return configure_args
