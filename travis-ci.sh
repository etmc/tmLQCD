#!/bin/bash

########################################################################
#
# Copyright (C) 2017 Martin Ueding <dev@martin-ueding.de>
#
# This file is part of tmLQCD.
#
# tmLQCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# tmLQCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

# Compiles and tests tmLQCD on the [Travis CI](https://travis-ci.org/) infrastructure.

set -e
set -u
set -x

# Compile C-LIME.
pushd ..
git clone https://github.com/usqcd-software/c-lime.git
pushd c-lime
./autogen.sh
./configure
make -j $(nproc)
popd
popd

if ! [[ -f git_hash.h ]]; then
    echo "#ifndef _GIT_HASH_H" > git_hash.h
    echo "#define _GIT_HASH_H" >> git_hash.h
    echo "const char git_hash[] = {\"travisbuild\"};" >> git_hash.h
    echo "#endif /* _GIT_HASH_H */" >> git_hash.h
fi

# Compile tmLQCD.
sudo apt-get update
sudo apt-get install -y flex libblas-dev liblapack-dev gfortran

autoconf

./configure \
    --disable-mpi \
    --with-lapack='-llapack -lblas' \
    --with-limedir=$PWD/../c-lime \
    CC=/usr/bin/gcc \
    CXX=/usr/bin/g++ \
    CFLAGS='-O2 --std=c99 -fopenmp -g -fPIC' \
    CXXFLAGS='-O2 --std=c++11 -fopenmp -g -fPIC' \
    LIBS='-fopenmp' \
    || ( echo; echo '###############################################################################'; echo '#                                Configure Log                                #'; echo '###############################################################################'; echo; set -x; cat config.log; exit 1)

make -j $(nproc)

# Run some tests.
cp sample-input/sample-hmc0.input travis.input
sed -i 's/Measurements = 1000/Measurements = 1/' travis.input
./hmc_tm -f travis.input
