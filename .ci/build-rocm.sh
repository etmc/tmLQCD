#!/bin/bash

set -xeuo pipefail

echo "VARIABLE = $VARIABLE"

export SPACK_SYSTEM_CONFIG_PATH=/user-environment/config
export CICD_SRC_DIR=$PWD
export QUDA_SRC_DIR=$PWD/deps/src/quda
export SPACK_PYTHON=$(which python3.6) # must be <=3.12, system python is 3.6

# QUDA git, branch and commit
export QUDA_GIT_REPO="${QUDA_GIT_REPO:=https://github.com/lattice/quda.git}"
export QUDA_GIT_BRANCH="${QUDA_GIT_BRANCH:=develop}"
export QUDA_GIT_COMMIT="${QUDA_GIT_COMMIT:=$(git ls-remote ${QUDA_GIT_REPO} refs/heads/${QUDA_GIT_BRANCH} | awk '{print $1}')}"

# make sure we keep the stage direcorty
spack config --scope=user add config:build_stage:/dev/shm/spack-stage
# we might need to install dependencies too, e.g. nlcglib in case of API changes
spack config --scope=user add config:install_tree:root:/dev/shm/spack-stage

spack env create -d ./spack-env

# add local repository with current tmlqcd recipe
spack -e ./spack-env repo add $REPO

spack -e ./spack-env config add "packages:all:variants:[amdgpu_target=${ROCM_ARCH},amdgpu_target_sram_ecc=${ROCM_ARCH},+rocm]"

spack -e ./spack-env add $SPEC

# for tmlqcd use local src instead of fetch git
spack -e ./spack-env develop -p ${CICD_SRC_DIR} tmlqcd@cicd

# for quda use local src instead of fetch git, to be able to tests against
# differnt repo, branch, commit and also to support that quda branch develop is
# a moving target
spack -e ./spack-env develop -p ${QUDA_SRC_DIR} quda@cicd

# display spack.yaml
cat ./spack-env/spack.yaml

spack -e ./spack-env concretize
spack -e ./spack-env install

# the tar pipe below expects a relative path
builddir=$(spack -e ./spack-env location -b tmlqcd)

# create a symlink to spack build directory (keep in artifacts)
tar -cf builddir.tar $builddir
