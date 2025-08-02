# MuFoLAB - Multi-Format Linear Algebra Benchmarks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14540573.svg)](https://doi.org/10.5281/zenodo.14540573)

This repository provides facilities for running automated large-scale
benchmarks on numerical linear algebra problems across multiple machine
number formats, including IEEE 754 floating-point numbers, bfloat16,
Posit arithmetic, and Takum arithmetic.

The repository includes several matrix test sets: sparse matrices are
sourced from the SuiteSparse Matrix Collection and stored in a compressed
format to facilitate deployment on HPC systems, while full matrices are
generated pseudo-randomly. The currently implemented algorithms encompass
both direct and indirect solvers, including UMFPACK-like LU, SPQR-like QR,
three-precision MPIR, and GMRES with iLU preconditioning.

There are also facilities to scrape graph matrices from the Network
Repository, yielding a test set of normalized, symmetrized graph
Laplacians. Based on this, interfaces for the computation of eigenvalues
and eigenvectors have been added.

## Setting up the environment

### Hardware Requirements

The machine must support 80-bit extended precision floating-point types
to properly simulate 64-bit tapered-precision arithmetic. This is the
case for x86_64, but not for ARM or RISC-V. At least 32 GB of RAM and
2 GB of storage (and an additional 1 GB if you are using the docker
shell environment) are required (including the compressed artifact
archive).

### Docker

This repository provides a recipe for a docker shell environment that
gets automatically set up with all necessary dependencies. Simply run

```sh
./run_docker_shell.sh
```

and you will be dropped in a shell of a docker container where the
repository is mounted as the working directory. Any files changed or
added in the container will persist after closing it, thus leaving
the shell is not an issue.

### Native

Using docker is optional and you can also run everything natively. Just
ensure that you have julia 1.11 or newer, TeX Live and latexmk installed.
Before running any of the commands below, open a julia REPL with the
`--project` flag within the repository folder and run the commands

```julia
using Pkg
Pkg.instantiate()
```

to install all the necessary dependencies. All subsequent julia calls
within MuFoLAB also use the `--project` flag; this way the installed
dependencies are bound to this repository and not installed globally.

## Getting started

You can automatically run all benchmarks via

```sh
make
```

and specifically

```sh
make eigen
```

to run the eigenvalue benchmarks (and generate the plot file
plots/eigen.pdf) and

```sh
make solve
```

to run the solver benchmarks respectively.

Runtime parameters (thread count, datasets) can be controlled by editing
`config.mk`. Noteworthy here is the `--reduced-test-data` flag that can
be uncommented and thus added to JULIA_SCRIPT_FLAGS. This reduces the
test data insofar that, for instance, a full run for the eigenvalue
benchmarks takes just an hour on a typical desktop computer. The full
runs can take multiple days to complete.

## Authors and License

MuFoLAB is developed by Laslo Hunhold and James Quinlan and licensed
under the ISC license. See LICENSE for copyright and license details.
