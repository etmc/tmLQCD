# Information about the CSCS CI/CD pipeline

This document describes the external pipeline executed through CSCS.

The pipeline can be triggered by commenting on a pull request with

```
cscs-ci run default  # runs the default pipeline
```

An automatic trigger on all merge-requests is currently disabled.

This pipeline has 2 stages: `build` and `test`.

The `build` stage builds a uenv image that includes all necessary compilers, MPI libraries and other dependecies to build QUDA and tmLQCD against QUDA. In this stage, QUDA is built correctly for the GH200 machine at CSCS with all required build flags for production runs. The uenv recipe can be found [here](.ci/uenv-recipes/tmlqcd/daint-gh200).

In the `test` stage, the aforementioned uenv image is loaded, tmLQCD is built and linked against the QUDA library that is inside the image. Finally a minimal HMC is executed and checked against some reference data.

## Force recompilation of quda

Remove the build cache:

```bash
/capstor/scratch/cscs/${USER}/uenv-cache/user-environment/build_cache/linux-sles15-neoverse_v2/gcc-13.2.0/quda-*
/capstor/scratch/cscs/${USER}/uenv-cache/user-environment/build_cache/linux-sles15-neoverse_v2-gcc-13.2.0-quda*
```

Or increment the the version counter tag in [.ci/include/cscs/00-variables.yml](.ci/include/cscs/00-variables.yml):

```yml
  UENV_TAG: <increment>
```

and commit.

## Useful links

* [Repo CSCS Pipline Entrypoint](.ci/cscs_default_pipeline.yml)
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
