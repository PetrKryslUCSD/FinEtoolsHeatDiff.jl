var documenterSearchIndex = {"docs":
[{"location":"guide/guide.html#Guide","page":"Guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"The FinEtools package is used here to solve heat conduction problems.","category":"page"},{"location":"guide/guide.html#Modules","page":"Guide","title":"Modules","text":"","category":"section"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"The package FinEtoolsHeatDiff has the following structure:","category":"page"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"Top-level:    FinEtoolsHeatDiff is the  top-level module.\nHeat conduction: AlgoHeatDiffModule (algorithms), FEMMHeatDiffModule, FEMMHeatDiffSurfModule  (FEM machines  to evaluate  the  matrix and vector quantities), MatHeatDiffModule  (heat diffusion material)","category":"page"},{"location":"guide/guide.html#Heat-conduction-FEM-machines","page":"Guide","title":"Heat  conduction FEM machines","text":"","category":"section"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"There is one for  the interior integrals  and one for  boundary integrals. The  machine for the interior integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"Evaluate the conductivity matrix.\nEvaluate the load vector corresponding to prescribed temperature.","category":"page"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"The machine for the boundary integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"Compute surface heat transfer  matrix.\nCompute  the heat load vector for surface heat transfer.\nCompute the heat load vector  corresponding to prescribed temperatures on the boundary  with surface heat transfer. moved up","category":"page"},{"location":"guide/guide.html#Algorithms","page":"Guide","title":"Algorithms","text":"","category":"section"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.","category":"page"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"Algorithms typically (not always) accept a single argument, modeldata, a dictionary of data, keyed by Strings. Algorithms  also return modeldata,  typically  including additional key/value pairs that represent the data computed by the algorithm.","category":"page"},{"location":"guide/guide.html#Heat-diffusion-algorithms","page":"Guide","title":"Heat diffusion algorithms","text":"","category":"section"},{"location":"guide/guide.html","page":"Guide","title":"Guide","text":"There is an implementation of an algorithm for steady-state heat conduction.","category":"page"},{"location":"index.html#FinEtools-(Finite-Element-tools)-Documentation","page":"Home","title":"FinEtools (Finite Element tools) Documentation","text":"","category":"section"},{"location":"index.html#Conceptual-guide","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"index.html#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"man/types.md\",\n    \"man/functions.md\",\n]\nDepth = 2","category":"page"},{"location":"man/types.html#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"man/types.html#FEM-machines","page":"Types","title":"FEM machines","text":"","category":"section"},{"location":"man/types.html#Heat-diffusion","page":"Types","title":"Heat diffusion","text":"","category":"section"},{"location":"man/types.html","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtoolsHeatDiff.FEMMHeatDiffModule, FinEtoolsHeatDiff.FEMMHeatDiffSurfModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsHeatDiff.FEMMHeatDiffModule.FEMMHeatDiff","page":"Types","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.FEMMHeatDiff","text":"FEMMHeatDiff{S<:AbstractFESet, F<:Function, M<:MatHeatDiff} <: AbstractFEMM\n\nType for heat diffusion finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsHeatDiff.FEMMHeatDiffModule.FEMMHeatDiff-Union{Tuple{M}, Tuple{F}, Tuple{S}, Tuple{IntegDomain{S, F}, M}} where {S<:AbstractFESet, F<:Function, M<:MatHeatDiff}","page":"Types","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.FEMMHeatDiff","text":"FEMMHeatDiff(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M<:MatHeatDiff}\n\nConstruct with the default orientation matrix (identity).\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.FEMMHeatDiffSurf","page":"Types","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.FEMMHeatDiffSurf","text":"FEMMHeatDiffSurf{S<:AbstractFESet, F<:Function} <: AbstractFEMM\n\nType for heat diffusion finite element modeling machine for boundary integrals.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#Material-models","page":"Types","title":"Material models","text":"","category":"section"},{"location":"man/types.html#Material-models-for-heat-diffusion","page":"Types","title":"Material models for heat diffusion","text":"","category":"section"},{"location":"man/types.html","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtools.MatModule, FinEtoolsHeatDiff.MatHeatDiffModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtools.MatModule.AbstractMat","page":"Types","title":"FinEtools.MatModule.AbstractMat","text":"AbstractMat\n\nAbstract type of material.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff","page":"Types","title":"FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff","text":"MatHeatDiff{MTAN<:Function, MUPD<:Function} <: AbstractMat\n\nType of material model for heat diffusion.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff-Tuple{Any, Any}","page":"Types","title":"FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff","text":"MatHeatDiff(thermal_conductivity, specific_heat)\n\nConstruct material model for heat diffusion.\n\nSupply the matrix of thermal conductivity constants.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff-Tuple{Any}","page":"Types","title":"FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff","text":"MatHeatDiff(thermal_conductivity)\n\nConstruct material model for heat diffusion.\n\nSupply the matrix of thermal conductivity constants.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions.html#FEM-machines","page":"Functions","title":"FEM machines","text":"","category":"section"},{"location":"man/functions.html#Heat-diffusion","page":"Functions","title":"Heat diffusion","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsHeatDiff.FEMMHeatDiffModule, FinEtoolsHeatDiff.FEMMHeatDiffSurfModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtools.FEMMBaseModule.inspectintegpoints-Union{Tuple{F}, Tuple{T}, Tuple{FEMMHeatDiff, NodalField{Float64}, NodalField{T}, NodalField{Float64}, Vector{Int64}, F, Any}, Tuple{FEMMHeatDiff, NodalField{Float64}, NodalField{T}, NodalField{Float64}, Vector{Int64}, F, Any, Any}} where {T<:Number, F<:Function}","page":"Functions","title":"FinEtools.FEMMBaseModule.inspectintegpoints","text":"inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{FFlt}, u::NodalField{T}, temp::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:heatflux; context...) where {T<:Number, F<:Function}\n\nInspect integration point quantities.\n\nArguments\n\ngeom - reference geometry field\nu - displacement field (ignored)\ntemp - temperature field\nfelist - indexes of the finite elements that are to be inspected:   The fes to be included are: fes[felist].\ncontext    - structure: see the update!() method of the material.\ninspector - function with the signature       idat = inspector(idat, j, conn, x, out, loc);  where   idat - a structure or an array that the inspector may          use to maintain some state,  for instance heat flux, j is the         element number, conn is the element connectivity, out is the         output of the update!() method,  loc is the location of the         integration point in the reference configuration.\n\nOutput\n\nThe updated inspector data is returned.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffModule.capacity-Union{Tuple{A}, Tuple{FEMMHeatDiff, A, NodalField{Float64}, NodalField{Float64}}} where A<:AbstractSysmatAssembler","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.capacity","text":"capacity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}\n\nCompute the capacity matrix.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffModule.conductivity-Union{Tuple{A}, Tuple{FEMMHeatDiff, A, NodalField{Float64}, NodalField{Float64}}} where A<:AbstractSysmatAssembler","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.conductivity","text":"conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}\n\nCompute the conductivity matrix.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffModule.energy-Tuple{FEMMHeatDiff, NodalField{Float64}, NodalField{Float64}}","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.energy","text":"energy(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt})\n\nCompute the \"energy\" integral over the interior domain.\n\nThe \"energy\" density is the dot product of the gradient of temperature and the heat flux.\n\nArguments\n\nself = model machine,\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffModule.nzebcloadsconductivity-Union{Tuple{A}, Tuple{FEMMHeatDiff, A, NodalField{Float64}, NodalField{Float64}}} where A<:AbstractSysvecAssembler","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.nzebcloadsconductivity","text":"nzebcloadsconductivity(self::FEMMHeatDiff, assembler::A,  geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}\n\nCompute load vector for nonzero EBC of prescribed temperature.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.nzebcsurfacetransferloads-Union{Tuple{A}, Tuple{FEMMHeatDiffSurf, A, NodalField{Float64}, NodalField{Float64}}} where A<:AbstractSysvecAssembler","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.nzebcsurfacetransferloads","text":"nzebcsurfacetransferloads(self::FEMMHeatDiffSurf, assembler::A,\n  geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}\n\nCompute load vector for nonzero EBC for fixed temperature.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field; the fixed_values attribute needs to list     the prescribed values of the temperature. This is performed with     applyebc!.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransfer-Union{Tuple{A}, Tuple{FEMMHeatDiffSurf, A, NodalField{Float64}, NodalField{Float64}}} where A<:AbstractSysmatAssembler","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransfer","text":"surfacetransfer(self::FEMMHeatDiffSurf,  assembler::A,\n  geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}\n\nCompute the surface heat transfer matrix.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransferloads-Union{Tuple{A}, Tuple{FEMMHeatDiffSurf, A, NodalField{Float64}, NodalField{Float64}, NodalField{Float64}}} where A<:AbstractSysvecAssembler","page":"Functions","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransferloads","text":"surfacetransferloads(self::FEMMHeatDiffSurf,  assembler::A,\n  geom::NodalField{FFlt}, temp::NodalField{FFlt},\n  ambtemp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}\n\nCompute the load vector corresponding to surface heat transfer.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\nambtemp = ambient temperature field on the surface\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Algorithms","page":"Functions","title":"Algorithms","text":"","category":"section"},{"location":"man/functions.html#Heat-conduction","page":"Functions","title":"Heat conduction","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsHeatDiff.AlgoHeatDiffModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsHeatDiff.AlgoHeatDiffModule.steadystate-Tuple{Dict{String, Any}}","page":"Functions","title":"FinEtoolsHeatDiff.AlgoHeatDiffModule.steadystate","text":"steadystate(modeldata::FDataDict)\n\nSteady-state heat conduction solver.\n\nArgument\n\nmodeldata = dictionary with items\n\n\"fens\"  = finite element node set\n\"regions\"  = array of region dictionaries\n\"essential_bcs\" = array of essential boundary condition dictionaries\n\"convection_bcs\" = array of convection boundary condition dictionaries\n\"flux_bcs\" = array of flux boundary condition dictionaries\n\nFor each region (connected piece of the domain made of a particular material), mandatory, the  region dictionary  contains items:\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"Q\" = material internal heat generation rate (optional; default  0.0)\n\nFor essential boundary conditions (optional) each dictionary would hold\n\n\"temperature\" = fixed (prescribed) temperature (scalar),  or         a function with signature             function T = f(x)         If not given, zero temperatures assumed.\n\"node_list\" = list of nodes on the boundary to which the condition applies         (mandatory)\n\nFor convection boundary conditions (optional) each dictionary may hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"ambient_temperature\" = fixed (prescribed) ambient temperature (scalar)     If not given, zero temperatures assumed.\n\nFor flux boundary conditions (optional) each dictionary would hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"normal_flux\" = normal component of the flux through the boundary (scalar)     Positive  when outgoing.\n\nOutput\n\nmodeldata= the dictionary on input is augmented with\n\n\"geom\" = the nodal field that is the geometry\n\"temp\" = the nodal field that is the computed temperature\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Material-models","page":"Functions","title":"Material models","text":"","category":"section"},{"location":"man/functions.html#Material-models-for-heat-diffusion","page":"Functions","title":"Material models for heat diffusion","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsHeatDiff.MatHeatDiffModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsHeatDiff.MatHeatDiffModule.tangentmoduli!","page":"Functions","title":"FinEtoolsHeatDiff.MatHeatDiffModule.tangentmoduli!","text":"tangentmoduli!(self::MatHeatDiff, kappabar::FFltMat, t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)\n\nCalculate the thermal conductivity matrix.\n\nkappabar = matrix of thermal conductivity (tangent moduli) in material coordinate system, supplied as a buffer and overwritten.\n\n\n\n\n\n","category":"function"},{"location":"man/functions.html#FinEtoolsHeatDiff.MatHeatDiffModule.update!","page":"Functions","title":"FinEtoolsHeatDiff.MatHeatDiffModule.update!","text":"update!(self::MatHeatDiff, heatflux::FFltVec, output::FFltVec, gradT::FFltVec, t::FFlt= 0.0, dt::FFlt= 0.0,  loc::FFltMat=FFltMat[], label::FInt=0, quantity=:nothing)\n\nUpdate material state.\n\nArguments\n\ngradT = thermal gradient vector,\nt = current time,\ndt = current time step,\nloc = location of the quadrature point in global Cartesian coordinates,\nlabel = label of the finite element in which the quadrature point is located.\nquantity = quantity to be output (:heatflux)\n\nOutput\n\nheatflux = heat flux vector, allocated by the caller with a size of the embedding space. The components of the heat flux vector are calculated and stored in the heatflux vector.\noutput =  array which is (if necessary) allocated  in an appropriate size, filled with the output quantity, and returned.\n\n\n\n\n\n","category":"function"}]
}