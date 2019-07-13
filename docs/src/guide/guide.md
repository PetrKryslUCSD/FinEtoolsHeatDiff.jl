# Guide

The [`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest/index.html) package is used here to solve heat conduction problems.

## Modules

The package `FinEtoolsHeatDiff` has the following structure:

- Top-level:
     `FinEtoolsHeatDiff` is the  top-level module.

- Heat conduction: `AlgoHeatDiffModule` (algorithms), `FEMMHeatDiffModule`, `FEMMHeatDiffSurfModule`  (FEM machines  to evaluate  the  matrix and vector quantities), `MatHeatDiffModule`  (heat diffusion material)


###  Heat  conduction FEM machines

There is one for  the interior integrals  and one for  boundary integrals.
The  machine for the interior integrals can be used to compute:

- Evaluate the conductivity matrix.

- Evaluate the load vector corresponding to prescribed temperature.

The machine for the boundary integrals can be used to compute:

- Compute surface heat transfer  matrix.

- Compute  the heat load vector for surface heat transfer.

- Compute the heat load vector  corresponding to prescribed temperatures on the boundary  with surface heat transfer.

## Material and Material Orientation

The material response  is described in  material-point-attached coordinate system. These coordinate systems  are Cartesian, and the material coordinate system is typically chosen to make  the response particularly simple.  So for orthotropic or transversely isotropic materials the axes would be aligned with the axes of orthotropy.

The type `CSys` (module `CSysModule`) is the updater of the material coordinate system matrix. The object is equipped with a callback to store the current orientation matrix. For instance: the coordinate system for an orthotropic material wound around a cylinder could be described in the coordinate system `CSys(3, 3, updatecs!)`, where the callback is defined as

```julia
function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    csmatout[:, 2] = [0.0 0.0 1.0]
    csmatout[:, 3] = XYZ
    csmatout[3, 3] = 0.0
    csmatout[:, 3] = csmatout[:, 3]/norm(csmatout[:, 3])
    csmatout[:, 1] = cross(csmatout[:, 2], csmatout[:, 3])
end
```

## Algorithms

Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.

Algorithms typically (not always) accept a single argument, `modeldata`, a dictionary of data, keyed by Strings. Algorithms  also return `modeldata`,  typically  including additional key/value pairs that represent the data computed by the algorithm.

### Heat diffusion algorithms

There is an implementation of an algorithm for steady-state heat conduction.
