"""
FinEtools (C) 2017-2023, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics. Package for heat diffusion problems.
"""
module FinEtoolsHeatDiff

__precompile__(true)

include("MatHeatDiffModule.jl")
include("FEMMHeatDiffModule.jl")
include("FEMMHeatDiffSurfModule.jl")
include("AlgoHeatDiffModule.jl")

# Exports follow:

###########################################################################
# Heat diffusion functionality
###########################################################################
using FinEtoolsHeatDiff.MatHeatDiffModule: MatHeatDiff
# Exported: type of heat-diffusion  material
export MatHeatDiff

using FinEtoolsHeatDiff.FEMMHeatDiffModule: FEMMHeatDiff,
    conductivity, energy, inspectintegpoints, capacity
# Exported: type  for linear heat diffusion and discretization methods
export FEMMHeatDiff, conductivity, energy, inspectintegpoints, capacity

using FinEtoolsHeatDiff.FEMMHeatDiffSurfModule: FEMMHeatDiffSurf,
    surfacetransfer, surfacetransferloads
# Exported: type  for linear heat diffusion boundary conditions and discretization methods
export FEMMHeatDiffSurf, surfacetransfer, surfacetransferloads

end # module
