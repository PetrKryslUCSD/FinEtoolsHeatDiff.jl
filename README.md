[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# FinEtoolsHeatDiff: Linear stress analysis application

`FinEtools` is a package for basic operations on finite element meshes.
`FinEtoolsHeatDiff` is a package using `FinEtools` to solve linear heat conduction (diffusion) problems.

## News

- 06/11/2019: Applications are now separated  out from the `FinEtools` package.

[Past news](oldnews.md)

## How to run

The [FinEtools](https://github.com/PetrKryslUCSD/FinEtools.jl) package is
needed. The entire setup of `FinEtoolsHeatDiff` can be performed with
```julia
] activate .; instantiate
```

The package `FinEtoolsHeatDiff` can be tested as
```julia
] activate .; instantiate; test
```

There are a number of examples, which may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
