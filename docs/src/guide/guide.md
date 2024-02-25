# Guide

The [`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest/index.html) package is used here to solve heat conduction problems.

Tutorials  are provided in the form of Julia scripts and Markdown files in a dedicated folder: [`index of tutorials`](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl/blob/main/tutorials/index.md). 


## Modules

The package `FinEtoolsHeatDiff` has the following structure:

- **Top-level**:
     `FinEtoolsHeatDiff` is the  top-level module.

- **Heat conduction**: `AlgoHeatDiffModule` (algorithms), `FEMMHeatDiffModule`, `FEMMHeatDiffSurfModule`  (FEM machines  to evaluate  the  matrix and vector quantities), `MatHeatDiffModule`  (heat diffusion material)


##  Heat  conduction FEM machines

There is one for  the interior integrals  and one for the boundary
integrals. The  machine for the interior integrals can be used to
compute:

- Conductivity matrix.

- Load vector corresponding to prescribed temperature.

The machine for the boundary integrals can be used to compute:

- Surface heat transfer  matrix.

- Heat load vector for surface heat transfer.


## Algorithms

Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.

Algorithms typically (not always) accept a single argument, `modeldata`, a dictionary of data, keyed by Strings. Algorithms  also return `modeldata`,  typically  including additional key/value pairs that represent the data computed by the algorithm.

### Heat diffusion algorithms

There is an implementation of an algorithm for steady-state heat conduction.

## Examples

### Two-dimensional heat transfer with convection: convergence study.


Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
short edge the temperature is fixed at 100 °C, and on one long edge the
plate is perfectly insulated so that the heat flux is zero through that
edge. The other two edges are losing heat via convection to an ambient
temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
.°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
There is no internal generation of heat. Calculate the temperature 0.2 m
along the un-insulated long side, measured from the intersection with the
fixed temperature side. The reference result is 18.25 °C.

The reference temperature at the point A  is 18.25 °C according to the
NAFEMS publication (which cites the book Carslaw, H.S. and J.C. Jaeger,
Conduction of Heat in Solids. 1959: Oxford University Press).

The present  tutorial will investigate the reference temperature  and it
will attempt to  estimate the  limit value more precisely using a
sequence of meshes and Richardson's extrapolation.

We begin by defining a few geometric parameters.
```julia
using FinEtools
# Geometrical dimensions
Width = 0.6 * phun("M")
Height = 1.0 * phun("M")
HeightA = 0.2 * phun("M")
Thickness = 0.1 * phun("M")
tolerance = Width / 1000
```

And now we are ready to define the conductivity of the material.
```julia
using FinEtoolsHeatDiff
# Conductivity matrix
kappa = [52.0 0; 0 52.0] * phun("W/(M*K)") 
m = MatHeatDiff(kappa)
```

The surface heat transfer coefficient (film coefficient) is also defined,
```julia
# Surface heat transfer coefficient
h = 750 * phun("W/(M^2*K)")
```
as is the temperature along one edge.
```julia
# Prescribed temperature.
T1 = 100 * phun("K")
```

Now we are ready to round the simulation five times, for progressively
increasing refinement levels, and collect the computed value of the reference
quantity.

```julia
# Five progressively refined models will be created and solved. 
modeldata = nothing
resultsTempA = Float64[]
params = Float64[]
for nref in 2:6
    # The mesh is created from two rectangular blocks to begin with.
    fens, fes = T3blockx([0.0, Width], [0.0, HeightA])
    fens2, fes2 = T3blockx([0.0, Width], [HeightA, Height])
    # The meshes are then glued into a single entity.
    fens, newfes1, fes2 = mergemeshes(fens, fes, fens2, fes2, tolerance)
    fes = cat(newfes1, fes2)
    # Refine the mesh desired number of times.
    for ref in 1:nref
        fens, fes = T3refine(fens, fes)
    end
    # The boundary is extracted.
    bfes = meshboundary(fes)
    # The prescribed temperature is applied along edge 1 (the bottom
    # edge in Figure 1).
    list1 = selectnode(fens; box = [0.0 Width 0.0 0.0], inflate = tolerance)
    essential1 = FDataDict("node_list" => list1, "temperature" => T1)
    # The convection (surface heat transfer) boundary condition is applied
    # along the edges 2,3,4. 
    list2 = selectelem(fens, bfes; box = [Width Width 0.0 Height], inflate = tolerance)
    list3 = selectelem(fens, bfes; box = [0.0 Width Height Height], inflate = tolerance)
    # The boundary integrals are evaluated using a surface FEMM.
    cfemm = FEMMHeatDiffSurf(
        IntegDomain(subset(bfes, vcat(list2, list3)), GaussRule(1, 3), Thickness),
        h,
    )
    convection1 = FDataDict("femm" => cfemm, "ambient_temperature" => 0.0)
    # The interior integrals are evaluated using a volume FEMM.
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(3), Thickness), m)
    region1 = FDataDict("femm" => femm)
    # Make the model data
    modeldata = FDataDict(
        "fens" => fens,
        "regions" => [region1],
        "essential_bcs" => [essential1],
        "convection_bcs" => [convection1],
    )
    # Call the solver
    modeldata = AlgoHeatDiffModule.steadystate(modeldata)
    # Locate the node at the point A  [coordinates (Width,HeightA)].
    list4 = selectnode(fens; box=[Width Width HeightA HeightA], inflate=tolerance)
    # Collect the temperature  at the point A.
    Temp = modeldata["temp"]
    println("$(Temp.values[list4][1])")
    push!(resultsTempA, Temp.values[list4][1])
    push!(params, 1.0 / 2^nref)
end
```

The computed results can be processed with Richardson extrapolation to arrive at an estimate of the true solution.
```julia
# These are the computed results for the temperature at point A:
println("$( resultsTempA  )")
# Richardson extrapolation can be used to estimate the limit.
solnestim, beta, c, residual = richextrapol(resultsTempA[(end-2):end], params[(end-2):end])
println("Solution estimate = $(solnestim)")
println("Convergence rate estimate  = $(beta)")
```

In order to visualize the results, we export to Paraview. The geometry is
two-dimensional: this means we can visualize the temperature as a three
dimensional surface raised above the mesh.

```julia
# Postprocessing
geom = modeldata["geom"]
Temp = modeldata["temp"]
regions = modeldata["regions"]
vtkexportmesh(
    "T4NAFEMS--T3-solution.vtk",
    connasarray(regions[1]["femm"].integdomain.fes),
    [geom.values (Temp.values / 100)],
    FinEtools.MeshExportModule.VTK.T3;
    scalars = [("Temperature", Temp.values)],
)
vtkexportmesh(
    "T4NAFEMS--T3-mesh.vtk",
    connasarray(regions[1]["femm"].integdomain.fes),
    geom.values,
    FinEtools.MeshExportModule.VTK.T3,
)
```

![Alt Visualization of the temperature field](T4NAFEMS--T3-solution.png)
