
Issues and ideas:

-- Documenter:
using FinEtoolsHeatDiff
using DocumenterTools
Travis.genkeys(user="PetrKryslUCSD", repo="https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl")

-- using JuliaFormatter
using JuliaFormatter
format("./src", SciMLStyle(), annotate_untyped_fields_with_any=true)   


-- FEMMHeatDiffSurf{ID <: IntegDomain, FT<:Number} <: AbstractFEMM
The surface heat transfer coefficient should probably be a data cache.
