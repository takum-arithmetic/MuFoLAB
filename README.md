# MuFoLAB - Multi-Format Linear Algebra Benchmarks

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

## Getting started

You can automatically run the benchmarks via

```sh
make
```

Runtime parameters (thread count, datasets) can be controlled by editing
`config.mk`.

## Authors and License

MuFoLAB is developed by Laslo Hunhold and James Quinlan and licensed
under the ISC license. See LICENSE for copyright and license details.