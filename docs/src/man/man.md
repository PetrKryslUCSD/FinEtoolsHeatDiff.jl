# Manual


```@meta
CurrentModule = FinEtoolsHeatDiff
```

## FEM machines

### Heat diffusion: volume

```@docs
FEMMHeatDiff
capacity
conductivity
energy
inspectintegpoints
```

### Heat diffusion: surface

```@docs
FEMMHeatDiffSurf
surfacetransfer
surfacetransferloads
```

## Algorithms

### Heat conduction

```@docs
AlgoHeatDiffModule.steadystate
```

## Material models

### Material models for heat diffusion

```@docs
MatHeatDiff
MatHeatDiffModule.tangentmoduli!
MatHeatDiffModule.update!
```

## Modules

```@docs
FinEtoolsHeatDiff.FinEtoolsHeatDiff
FinEtoolsHeatDiff.FEMMHeatDiffModule
FinEtoolsHeatDiff.FEMMHeatDiffSurfModule
FinEtoolsHeatDiff.AlgoHeatDiffModule
FinEtoolsHeatDiff.MatHeatDiffModule
```

