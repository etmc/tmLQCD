# Information about the CSCS CI/CD pipeline

This document describes the external pipeline executed through CSCS.

The pipeline can be triggered by commenting on a pull request with

```
cscs-ci run default  # runs the default pipeline (on GH200 nodes @ CSCS)
cscs-ci run beverin  # runs the beverin pipeline (on MI300A nodes @ CSCS)
```

An automatic trigger on all merge-requests is currently disabled.

This pipeline has 3 stages: `prepare`, `build` and `test`.

## `prepare` stage

The `prepare` stage builds an uenv image that includes all necessary compilers, MPI libraries and other dependecies to build QUDA and tmLQCD against QUDA. The uenv recipe can be found [here for GH200](uenv-recipes/tmlqcd/daint-gh200) and [here for MI300A](uenv-recipes/tmlqcd/beverin-mi300).

## `build` stage

In the `build` stage, the aforementioned uenv image is loaded, tmLQCD and QUDA are built using their spack packages using the dependencies from the base image. This stage exposes an artifact with tmLQCD/QUDA binaries. For tmLQCD, the current branch is compiled. For QUDA the following environment variables are respected:

  * `QUDA_GIT_REPO`: the git repository URL to use as source (defaults to `https://github.com/lattice/quda.git`)
  * `QUDA_GIT_BRANCH`: the git branch to compile (defaults to `develop`)
  * `QUDA_GIT_COMMIT`: the git commit to compile (defaults to the current head commit of `QUDA_GIT_BRANCH`)

Then QUDA is cloned and compiled, completely bypassing the spack compile cache.

## `test` stage

In the `test` stage, the aforementioned uenv image is loaded, tmLQCD and QUDA are unpacked from the artifact. Finally a minimal HMC is executed and checked against some reference data.

## Force recompilation of base image in `prepare` stage

Remove the build cache:

```bash
/capstor/scratch/cscs/${USER}/uenv-cache/user-environment/build_cache/linux-sles15-neoverse_v2/gcc-13.2.0/tmlqcd-*
/capstor/scratch/cscs/${USER}/uenv-cache/user-environment/build_cache/linux-sles15-neoverse_v2-gcc-13.2.0-tmlqcd*
```

Or increment the the version counter tag in [.ci/include/cscs/00-variables.yml](include/cscs/00-variables.yml):

```yml
  UENV_TAG: <increment>
```

and commit.

## Useful links

* [Repo CSCS Pipline Entrypoint](cscs_default_pipeline.yml)
* [Repo Webhooks](https://github.com/etmc/tmLQCD/settings/hooks)
* [CSCS Knowledge Base](https://confluence.cscs.ch/display/KB)
* [CSCS Knowledge Base / CI/CD](https://confluence.cscs.ch/pages/viewpage.action?pageId=868812112)
* [CSCS Registered Projects](https://cicd-ext-mw.cscs.ch/ci/overview)
* [CSCS Developer Portal](https://developer.cscs.ch/)
* [CSCS Alps Uenv Recipes](https://github.com/eth-cscs/alps-uenv)
* [CSCS Uenv Documentation](https://eth-cscs.github.io/uenv/)
* [CSCS Uenv Writing Documentation](https://eth-cscs.github.io/alps-uenv/)
* [CSCS Status Page](https://status.cscs.ch/)
* [CSCS Spack Base Containers](https://github.com/orgs/eth-cscs/packages/container/package/docker-ci-ext%2Fspack-base-containers%2Fspack-build)
* [Sirius CI/CD](https://github.com/electronic-structure/SIRIUS/tree/develop/ci) where this one is based upon
