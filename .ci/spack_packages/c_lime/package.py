# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems import cmake
from spack_repo.builtin.build_systems.cmake import CMakePackage, generator

from spack.package import *


class CLime(CMakePackage):
    """LIME (which can stand for Lattice QCD Interchange Message Encapsulation
    or more generally, Large Internet Message Encapsulation) is a simple
    packaging scheme for combining records containing ASCII and/or binary
    data."""

    homepage = "https://usqcd-software.github.io/c-lime/"
    url = "https://github.com/usqcd-software/c-lime/archive/qio2-3-9.tar.gz"
    git = "https://github.com/usqcd-software/c-lime.git"

    license("GPL-2.0-or-later")

    version(
        "2-3-9",
        sha256="7b9aeadd4dfec50e24da3e7e729f56abf95c9192612c41515fe27b2158773aac",
    )
    version("master", branch="master")

    variant("pic", default=True, description="Enable position independent code")
    variant("shared", default=True, description="Enable shared libraries")

    depends_on("c", type="build")  # generated


class CMakeBuilder(cmake.CMakeBuilder):
    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
    return args
