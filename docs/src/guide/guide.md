# Guide

The [`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest/index.html) package is used here to solve heat conduction problems.

## Modules

The package `FinEtoolsHeatDiff` has the following structure:

- **Top-level**:
     `FinEtoolsHeatDiff` is the  top-level module.

- **Heat conduction**: `AlgoHeatDiffModule` (algorithms), `FEMMHeatDiffModule`, `FEMMHeatDiffSurfModule`  (FEM machines  to evaluate  the  matrix and vector quantities), `MatHeatDiffModule`  (heat diffusion material)


##  Heat  conduction FEM machines

There is one for  the interior integrals  and one for  boundary
integrals. The  machine for the interior integrals can be used to
compute:

- Evaluate the conductivity matrix.

- Evaluate the load vector corresponding to prescribed temperature.

The machine for the boundary integrals can be used to compute:

- Compute surface heat transfer  matrix.

- Compute  the heat load vector for surface heat transfer.

- Compute the heat load vector  corresponding to prescribed temperatures on the
  boundary  with surface heat transfer. moved up

## Algorithms

Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.

Algorithms typically (not always) accept a single argument, `modeldata`, a dictionary of data, keyed by Strings. Algorithms  also return `modeldata`,  typically  including additional key/value pairs that represent the data computed by the algorithm.

### Heat diffusion algorithms

There is an implementation of an algorithm for steady-state heat conduction.


