The software ships with a CMake environment, which will configure and build the
programmes. It is recommended to configure and build the executables in a
separate build directory. This also allows to have several builds with different
options from the same source code directory.

## Prerequisites

In order to compile the programmes the `LAPACK` library (fortran version) needs to be installed. CMake will search for the
library in all default directories. Also the latest version (tested is version
1.2.3) of `C-LIME` must be available, which is used as
a packaging scheme to read and write gauge configurations and propagators to
files.

## Configuring the hmc package
:label{sec:config}

The build system uses CMake to configure and build the hmc package. The
following list gives all options (OFF by default unless specified):
- `CMAKE_POSITION_INDEPENDENT_CODE`: Build a position independent
  code. **ON** by default.
- `BUILD_SHARED_LIBS`: Build the shared version of the hmc library.
- `TM_USE_FFTW`: Enable fftw support. 
- `TM_USE_CUDA`: Enable CUDA support.
- `TM_USE_HIP`: Enable HIP support (AMD or NVidia GPUs)
- `TM_USE_DDalphaAMG`: Enable DDalphaAMG support.
- `TM_USE_LEMON`: Use the lemon io library.
- `TM_USE_OMP`: Enable OpenMP (**ON** by default)
- `TM_FIXEDVOLUME`: Fix volume at compile time.
- `TM_ENABLE_ALIGNMENT`: Automatically or expliclty align arrays to
  byte number. auto, none, 16, 32, 64.
- `TM_USE_GAUGE_COPY`: Enable use of a copy of the gauge field (**ON**
  by default). See section ref{sec:dirac} for details on this option. It will
  increase the memory requirement of the code.
- `TM_USE_HALFSPINOR`: Use a Dirac Op. with halfspinor exchange (**ON**
  by default). See sub-section ref{sec:dirac} for details. 
- `TM_USE_QUDA`: Enable QUDA support.
- `TM_USE_SHMEM`: Use shmem API.
- `TM_ENABLE_WARNINGS`: Enable all warnings (**ON** by default).
- `TM_ENABLE_TESTS`: Enable tests.
- `TM_USE_QPHIX`: Enable QPhiX.
  - `TM_QPHIX_SOALEN`: QPhiX specific parameter (default is 4)
  - **QPHIX_DIR**: Directory where QPhiX is installed.
    The QPhiX current CMake build system does not export all information (
    include and lib directories) that are needed to compile hmc.
  - **QMP_DIR**: Directory where QMP is installed (
    QPhiX dependency).
    The QPhiX current CMake build system does not export all information about the
    include and lib directories nor its dependencies (QMP in that case).
- `TM_USE_MPI`: Enable MPI support.
  - `TM_PERSISTENT_MPI`: Use persistent MPI calls for halfspinor.
  - `TM_NONBLOCKING_MPI`: Use non-blocking MPI calls for spinor and
    gauge.
  - `TM_MPI_DIMENSION`: Use $n$ dimensional parallelisation ($XYZT$)
    [default=4]. The number of parallel directions can be specified. $1, 2, 3$ and $4$
    dimensional parallelisation is supported.
  - `TM_USE_LEMON` Use the lemon io library

The following minimal list of commands will configure and build the hmc package with
minimal dependencies

```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/my_path -DCMAKE_PREFIX_PATH=/my_c_line_path ..
make -j
make install
```

These instructions assume that the `c-lime` package is installed in `/my_c_line_path`. By default `CMAKE_PREFIX_PATH` variable is a list
of paths separated by a semi-colunm containing the path of all installed to
dependencies.

Adding `-DTM_USE_MPI=ON` will enable MPI support with parallelization
over spatial and temporal dimensions. The command line is then

```bash
cmake -DCMAKE_INSTALL_PREFIX=/my_path -DCMAKE_PREFIX_PATH=/my_c_line_path -DTM_USE_MPI=ON ..
```

We can combine it with the lemon-io library (isntalled in `/my_lemon_path`)

```bash
cmake -DCMAKE_INSTALL_PREFIX=/my_path \
      -DCMAKE_PREFIX_PATH="/my_c_line_path;/my_lemon_path" \
      -DTM_USE_MPI=ON \
      -DTM_USE_LEMON=ON ..
```

`QUDA` support (installed in `/my_quda_path`) can be added with

```bash
cmake -DCMAKE_INSTALL_PREFIX=/my_path \
      -DCMAKE_PREFIX_PATH="/my_c_line_path;/my_lemon_path;/my_quda_path" \
      -DTM_USE_MPI=ON \
      -DTM_USE_LEMON=ON \
      -DTM_USE_QUDA \
      -DTM_USE_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=90 ..
```

Note that the command assumes that QUDA is compiled with `CUDA` support. AMD GPU
are also supported after replacing `-DTM_USE_CUDA=ON` with
`-DTM_USE_HIP=ON` and compiling `QUDA` with `HIP` support. The ROCM architecture is defined by the variable
`CMAKE_HIP_ARCHITECTURES=gfxxxx`.  An extra parameter `-DCMAKE_CXX_COMPILER=clang++` is needed because `QUDA` use the `ROCM clang++` 
compiler internally and the build will fail if `gcc` or any other compiler is used during 
link time. This option only affects the linking behavior not the compilation. The cmake command line for HIP/ROCM support is then
```bash
cmake -DCMAKE_INSTALL_PREFIX=/my_path \
    -DCMAKE_PREFIX_PATH="/my_c_line_path;/my_lemon_path;/my_quda_path" \
    -DTM_USE_MPI=ON \
    -DTM_USE_LEMON=ON \
    -DTM_USE_QUDA \
    -DTM_USE_HIP=ON \
    -DCMAKE_HIP_ARCHITECTURES=gfx90a \
    -DCMAKE_CXX_COMPILER=/opr/rocm/bin/clang++ ..
```

`QPhiX` and/or `DDalphaAMG` support can be added with

```bash
cmake -DCMAKE_INSTALL_PREFIX=/my_path \
      -DCMAKE_PREFIX_PATH="/my_c_line_path;/my_lemon_path;/my_quda_path;/my_path_ddalphaamg" \
      -DTM_USE_MPI=ON \
      -DTM_USE_LEMON=ON \
      -DTM_USE_QUDA=ON \
      -DTM_USE_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=90 \
      -DTM_USE_QPHIX=ON \
      -DQPHIX_DIR=/my_qphix_dir \
      -DTM_USE_DDalphaAMG=ON \
      -DQMP_DIR=/my_qmp_dir \
      -DTM_USE_OMP=ON ..
```

`QPhiX` cmake config support is incomplete and requires both the `QPhiX`
and `QMP` installation directories to work properly.

`CMake` has several relevant specific options that control the build. Compiler
options are defined by the variable `CMAKE_C_FLAGS` and `CMAKE_CXX_FLAGS`. CUDA and HIP compilations options are controlled by their
equivalent `CMAKE_{CUDA/HIP}_FLAGS`.

Adding for instance `-GNinja` to the `CMake` command line will use
ninja instead of make.
