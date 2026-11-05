#!/bin/bash

set -xeuo pipefail

# Shared, per-pipeline work dir on a filesystem visible to all jobs.
# Replaces builddir.tar: later jobs recompute the same path from CI_PIPELINE_ID.
# CI_PIPELINE_ID is instance-wide unique (CI_PIPELINE_IID is per project only).
export BUILD_DIR="${BUILD_DIR:-/capstor/store/cscs/pasc/c36/ci-stack/$CI_PIPELINE_ID}"
mkdir -p "$BUILD_DIR"

export SPACK_SYSTEM_CONFIG_PATH="/user-environment/config"
export SPACK_PYTHON="$(which python3.6)" # must be <=3.12, system python is 3.6

# dev-spec sources live on shared storage so their spack-build-* dirs persist
export CICD_SRC_DIR="$BUILD_DIR/tmlqcd"
export QUDA_SRC_DIR="$BUILD_DIR/quda"
export SPACK_ENV_DIR="$BUILD_DIR/spack-env"

# QUDA git, branch and commit
export QUDA_GIT_REPO="${QUDA_GIT_REPO:=https://github.com/lattice/quda.git}"
export QUDA_GIT_BRANCH="${QUDA_GIT_BRANCH:=develop}"
export QUDA_GIT_COMMIT="${QUDA_GIT_COMMIT:=$(git ls-remote ${QUDA_GIT_REPO} refs/heads/${QUDA_GIT_BRANCH} | awk '{print $1}')}"

# tmlqcd: copy the checkout to shared storage (retries reuse the same dir)
rm -rf "$CICD_SRC_DIR"
mkdir -p "$CICD_SRC_DIR"
cp -r "$PWD/." "$CICD_SRC_DIR/"   # or: rsync -a --exclude=.git "$PWD/" "$CICD_SRC_DIR/"

# obtain QUDA
rm -rf "$QUDA_SRC_DIR"
git clone -b "${QUDA_GIT_BRANCH}" "${QUDA_GIT_REPO}" "${QUDA_SRC_DIR}"
git -C "${QUDA_SRC_DIR}" checkout "${QUDA_GIT_COMMIT}"

rm -rf "$SPACK_ENV_DIR"
spack env create -d "$SPACK_ENV_DIR"

# Config goes into the environment (not --scope=spack, which would only live in
# the per-job spack clone), so later jobs using this env see the same paths.
# Stage is throwaway; the install tree must persist.
spack -e "$SPACK_ENV_DIR" config add "config:build_stage:/dev/shm/spack-stage"
spack -e "$SPACK_ENV_DIR" config add "config:install_tree:root:$BUILD_DIR/install"

# add local repository with current tmlqcd recipe
spack -e "$SPACK_ENV_DIR" repo add "${REPO}"

spack -e "$SPACK_ENV_DIR" config add "packages:all:variants:[${VARIANTS}]"

spack -e "$SPACK_ENV_DIR" add "${SPEC}"

# for tmlqcd use local src instead of fetch git
spack -e "$SPACK_ENV_DIR" develop -p "${CICD_SRC_DIR}" tmlqcd@cicd

# for quda use local src instead of fetch git, to be able to test against
# different repo, branch, commit and also to support that quda branch develop is
# a moving target
spack -e "$SPACK_ENV_DIR" develop -p "${QUDA_SRC_DIR}" quda@cicd

# display spack.yaml
cat "$SPACK_ENV_DIR/spack.yaml"

spack -e "$SPACK_ENV_DIR" concretize
spack -e "$SPACK_ENV_DIR" install --keep-prefix

# No tar step: the build dirs stay on shared storage at
#   $(spack -e "$SPACK_ENV_DIR" location -b tmlqcd)
#   $(spack -e "$SPACK_ENV_DIR" location -b quda)
# Downstream jobs use $BUILD_DIR directly.
spack -e "$SPACK_ENV_DIR" location -b tmlqcd
spack -e "$SPACK_ENV_DIR" location -b quda
