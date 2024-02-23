var documenterSearchIndex = {"docs":
[{"location":"guide/guide/#Guide","page":"Guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"The FinEtools package is used here to solve heat conduction problems.","category":"page"},{"location":"guide/guide/#Modules","page":"Guide","title":"Modules","text":"","category":"section"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"The package FinEtoolsHeatDiff has the following structure:","category":"page"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"Top-level:    FinEtoolsHeatDiff is the  top-level module.\nHeat conduction: AlgoHeatDiffModule (algorithms), FEMMHeatDiffModule, FEMMHeatDiffSurfModule  (FEM machines  to evaluate  the  matrix and vector quantities), MatHeatDiffModule  (heat diffusion material)","category":"page"},{"location":"guide/guide/#Heat-conduction-FEM-machines","page":"Guide","title":"Heat  conduction FEM machines","text":"","category":"section"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"There is one for  the interior integrals  and one for the boundary integrals. The  machine for the interior integrals can be used to compute:","category":"page"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"Conductivity matrix.\nLoad vector corresponding to prescribed temperature.","category":"page"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"The machine for the boundary integrals can be used to compute:","category":"page"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"Surface heat transfer  matrix.\nHeat load vector for surface heat transfer.","category":"page"},{"location":"guide/guide/#Algorithms","page":"Guide","title":"Algorithms","text":"","category":"section"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.","category":"page"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"Algorithms typically (not always) accept a single argument, modeldata, a dictionary of data, keyed by Strings. Algorithms  also return modeldata,  typically  including additional key/value pairs that represent the data computed by the algorithm.","category":"page"},{"location":"guide/guide/#Heat-diffusion-algorithms","page":"Guide","title":"Heat diffusion algorithms","text":"","category":"section"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"There is an implementation of an algorithm for steady-state heat conduction.","category":"page"},{"location":"guide/guide/#Example","page":"Guide","title":"Example","text":"","category":"section"},{"location":"guide/guide/","page":"Guide","title":"Guide","text":"(Image: Alt Visualization of the temperature field)","category":"page"},{"location":"#FinEtools-(Finite-Element-tools)-Documentation","page":"Home","title":"FinEtools (Finite Element tools) Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Conceptual-guide","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"man/man.md\",\n]\nDepth = 2","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"man/man/#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"man/man/","page":"Manual","title":"Manual","text":"CurrentModule = FinEtoolsHeatDiff","category":"page"},{"location":"man/man/#FEM-machines","page":"Manual","title":"FEM machines","text":"","category":"section"},{"location":"man/man/#Heat-diffusion:-volume","page":"Manual","title":"Heat diffusion: volume","text":"","category":"section"},{"location":"man/man/","page":"Manual","title":"Manual","text":"FEMMHeatDiff\ncapacity\nconductivity\nenergy\ninspectintegpoints","category":"page"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffModule.FEMMHeatDiff","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.FEMMHeatDiff","text":"FEMMHeatDiff{ID<:IntegDomain, M<:MatHeatDiff} <: AbstractFEMM\n\nType for heat diffusion finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffModule.capacity","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.capacity","text":"capacity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{GFT},  temp::NodalField{FT}) where {A<:AbstractSysmatAssembler, GFT, FT}\n\nCompute the capacity matrix.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"function"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffModule.conductivity","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.conductivity","text":"conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{GFT},  temp::NodalField{FT}) where {A<:AbstractSysmatAssembler, GFT, FT}\n\nCompute the conductivity matrix.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"function"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffModule.energy","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule.energy","text":"energy(self::FEMMHeatDiff, geom::NodalField{GFT},  temp::NodalField{FT}) where {GFT, FT}\n\nCompute the \"energy\" integral over the interior domain.\n\nThe \"energy\" density is the dot product of the gradient of temperature and the heat flux.\n\nArguments\n\nself = model machine,\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"function"},{"location":"man/man/#FinEtools.FEMMBaseModule.inspectintegpoints","page":"Manual","title":"FinEtools.FEMMBaseModule.inspectintegpoints","text":"inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{GFT}, u::NodalField{T}, temp::NodalField{FT}, felist::VecOrMat{IntT}, inspector::F, idat, quantity=:heatflux; context...) where {T<:Number, GFT, FT, IntT, F<:Function}\n\nInspect integration point quantities.\n\nArguments\n\ngeom - reference geometry field\nu - displacement field (ignored)\ntemp - temperature field\nfelist - indexes of the finite elements that are to be inspected:   The fes to be included are: fes[felist].\ncontext    - structure: see the update!() method of the material.\ninspector - function with the signature       idat = inspector(idat, j, conn, x, out, loc);  where   idat - a structure or an array that the inspector may          use to maintain some state,  for instance heat flux, j is the         element number, conn is the element connectivity, out is the         output of the update!() method,  loc is the location of the         integration point in the reference configuration.\n\nOutput\n\nThe updated inspector data is returned.\n\n\n\n\n\n","category":"function"},{"location":"man/man/#Heat-diffusion:-surface","page":"Manual","title":"Heat diffusion: surface","text":"","category":"section"},{"location":"man/man/","page":"Manual","title":"Manual","text":"FEMMHeatDiffSurf\nsurfacetransfer\nsurfacetransferloads","category":"page"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.FEMMHeatDiffSurf","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.FEMMHeatDiffSurf","text":"FEMMHeatDiffSurf{ID <: IntegDomain, FT<:Number} <: AbstractFEMM\n\nType for heat diffusion finite element modeling machine for boundary integrals.\n\n\n\n\n\n","category":"type"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransfer","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransfer","text":"surfacetransfer(self::FEMMHeatDiffSurf,  assembler::A, geom::NodalField{GFT}, temp::NodalField{FT})  where {A<:AbstractSysmatAssembler, GFT, FT}\n\nCompute the surface heat transfer matrix.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\n\n\n\n\n\n","category":"function"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransferloads","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule.surfacetransferloads","text":"surfacetransferloads(self::FEMMHeatDiffSurf,  assembler::A,  geom::NodalField{GFT}, temp::NodalField{FT},  ambtemp::NodalField{FT}) where {A<:AbstractSysvecAssembler, GFT, FT}\n\nCompute the load vector corresponding to surface heat transfer.\n\nArguments\n\nself = model machine,\nassembler = matrix assembler\ngeom = geometry field,\ntemp = temperature field\nambtemp = ambient temperature field on the surface\n\n\n\n\n\n","category":"function"},{"location":"man/man/#Algorithms","page":"Manual","title":"Algorithms","text":"","category":"section"},{"location":"man/man/#Heat-conduction","page":"Manual","title":"Heat conduction","text":"","category":"section"},{"location":"man/man/","page":"Manual","title":"Manual","text":"AlgoHeatDiffModule.steadystate","category":"page"},{"location":"man/man/#FinEtoolsHeatDiff.AlgoHeatDiffModule.steadystate","page":"Manual","title":"FinEtoolsHeatDiff.AlgoHeatDiffModule.steadystate","text":"steadystate(modeldata::FDataDict)\n\nSteady-state heat conduction solver.\n\nArgument\n\nmodeldata = dictionary with items\n\n\"fens\"  = finite element node set\n\"regions\"  = array of region dictionaries\n\"essential_bcs\" = array of essential boundary condition dictionaries\n\"convection_bcs\" = array of convection boundary condition dictionaries\n\"flux_bcs\" = array of flux boundary condition dictionaries\n\nFor each region (connected piece of the domain made of a particular material), mandatory, the  region dictionary  contains items:\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"Q\" = material internal heat generation rate (optional; default  0.0)\n\nFor essential boundary conditions (optional) each dictionary would hold\n\n\"temperature\" = fixed (prescribed) temperature (scalar),  or         a function with signature             function T = f(x)         If not given, zero temperatures assumed.\n\"node_list\" = list of nodes on the boundary to which the condition applies         (mandatory)\n\nFor convection boundary conditions (optional) each dictionary may hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"ambient_temperature\" = fixed (prescribed) ambient temperature (scalar)     If not given, zero temperatures assumed.\n\nFor flux boundary conditions (optional) each dictionary would hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"normal_flux\" = normal component of the flux through the boundary (scalar)     Positive  when outgoing.\n\nOutput\n\nmodeldata= the dictionary on input is augmented with\n\n\"geom\" = the nodal field that is the geometry\n\"temp\" = the nodal field that is the computed temperature\n\n\n\n\n\n","category":"function"},{"location":"man/man/#Material-models","page":"Manual","title":"Material models","text":"","category":"section"},{"location":"man/man/#Material-models-for-heat-diffusion","page":"Manual","title":"Material models for heat diffusion","text":"","category":"section"},{"location":"man/man/","page":"Manual","title":"Manual","text":"MatHeatDiff\nMatHeatDiffModule.tangentmoduli!\nMatHeatDiffModule.update!","category":"page"},{"location":"man/man/#FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff","page":"Manual","title":"FinEtoolsHeatDiff.MatHeatDiffModule.MatHeatDiff","text":"MatHeatDiff{FT, MTAN<:Function, MUPD<:Function} <: AbstractMat\n\nType of material model for heat diffusion.\n\n\n\n\n\n","category":"type"},{"location":"man/man/#FinEtoolsHeatDiff.MatHeatDiffModule.tangentmoduli!","page":"Manual","title":"FinEtoolsHeatDiff.MatHeatDiffModule.tangentmoduli!","text":"tangentmoduli!(self::MatHeatDiff, kappabar::Matrix{FT}, t = zero(FT), dt = zero(FT), loc::Matrix{FT} = reshape(FT[],0,0), label = 0) where {FT}\n\nCalculate the thermal conductivity matrix.\n\nkappabar = matrix of thermal conductivity (tangent moduli) in material coordinate system, supplied as a buffer and overwritten.\n\n\n\n\n\n","category":"function"},{"location":"man/man/#FinEtoolsHeatDiff.MatHeatDiffModule.update!","page":"Manual","title":"FinEtoolsHeatDiff.MatHeatDiffModule.update!","text":"update!(self::MatHeatDiff, heatflux::Vector{FT}, output::Vector{FT}, gradT::Vector{FT}, t= zero(FT), dt= zero(FT), loc::Matrix{FT}=reshape(FT[],0,0), label=0, quantity=:nothing) where {FT}\n\nUpdate material state.\n\nArguments\n\ngradT = thermal gradient vector,\nt = current time,\ndt = current time step,\nloc = location of the quadrature point in global Cartesian coordinates,\nlabel = label of the finite element in which the quadrature point is located.\nquantity = quantity to be output (:heatflux)\n\nOutput\n\nheatflux = heat flux vector, allocated by the caller with a size of the embedding space. The components of the heat flux vector are calculated and stored in the heatflux vector.\noutput =  array which is (if necessary) allocated  in an appropriate size, filled with the output quantity, and returned.\n\n\n\n\n\n","category":"function"},{"location":"man/man/#Modules","page":"Manual","title":"Modules","text":"","category":"section"},{"location":"man/man/","page":"Manual","title":"Manual","text":"FinEtoolsHeatDiff.FinEtoolsHeatDiff\nFinEtoolsHeatDiff.FEMMHeatDiffModule\nFinEtoolsHeatDiff.FEMMHeatDiffSurfModule\nFinEtoolsHeatDiff.AlgoHeatDiffModule\nFinEtoolsHeatDiff.MatHeatDiffModule","category":"page"},{"location":"man/man/#FinEtoolsHeatDiff.FinEtoolsHeatDiff","page":"Manual","title":"FinEtoolsHeatDiff.FinEtoolsHeatDiff","text":"FinEtools (C) 2017-2023, Petr Krysl\n\nFinite Element tools.  Julia implementation  of the finite element method for continuum mechanics. Package for heat diffusion problems.\n\n\n\n\n\n","category":"module"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffModule","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffModule","text":"FEMMHeatDiffModule\n\nModule for operations on interiors of domains to construct system matrices and system vectors for linear heat conduction/diffusion.\n\n\n\n\n\n","category":"module"},{"location":"man/man/#FinEtoolsHeatDiff.FEMMHeatDiffSurfModule","page":"Manual","title":"FinEtoolsHeatDiff.FEMMHeatDiffSurfModule","text":"FEMMHeatDiffSurfModule\n\nModule for operations on boundaries of domains to construct system matrices and system vectors for linear heat diffusion/conduction.\n\n\n\n\n\n","category":"module"},{"location":"man/man/#FinEtoolsHeatDiff.AlgoHeatDiffModule","page":"Manual","title":"FinEtoolsHeatDiff.AlgoHeatDiffModule","text":"AlgoHeatDiffModule\n\nModule for algorithms in linear heat conduction/diffusion  models.\n\n\n\n\n\n","category":"module"},{"location":"man/man/#FinEtoolsHeatDiff.MatHeatDiffModule","page":"Manual","title":"FinEtoolsHeatDiff.MatHeatDiffModule","text":"MatHeatDiffModule\n\nModule for linear heat diffusion material models.\n\n\n\n\n\n","category":"module"}]
}
