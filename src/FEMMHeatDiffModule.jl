"""
    FEMMHeatDiffModule

Module for operations on interiors of domains to construct system
matrices and system vectors for linear heat conduction/diffusion.
"""
module FEMMHeatDiffModule

using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, nodesperelem, manifdim, gradN!
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtools.CSysModule: CSys, updatecsmat!, csmat
using FinEtools.FieldModule: ndofs,
    gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
using FinEtools.NodalFieldModule: NodalField
using FinEtools.ElementalFieldModule: ElementalField
using FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!,
    SysvecAssembler
using FinEtools.ForceIntensityModule: ForceIntensity
using FinEtools.FEMMBaseModule: AbstractFEMM, finite_elements
import FinEtools.FEMMBaseModule: inspectintegpoints
using FinEtools.MatrixUtilityModule: add_gkgt_ut_only!,
    complete_lt!, locjac!, add_nnt_ut_only!, mulCAtB!, mulCAB!
using LinearAlgebra: norm, dot
using FinEtoolsHeatDiff.MatHeatDiffModule: MatHeatDiff, tangentmoduli!, update!
using FinEtools.FEMMBaseModule: bilform_diffusion, bilform_dot
using FinEtools.DataCacheModule: DataCache

"""
    FEMMHeatDiff{ID<:IntegDomain, M<:MatHeatDiff} <: AbstractFEMM

Type for heat diffusion finite element modeling machine.
"""
mutable struct FEMMHeatDiff{ID <: IntegDomain, CS <: CSys, M <: MatHeatDiff} <: AbstractFEMM
    integdomain::ID # geometry data
    mcsys::CS # updater of the material orientation matrix
    material::M # material object
end

"""
    FEMMHeatDiff(integdomain::ID, material::M) where {ID<:IntegDomain, M<:MatHeatDiff}

Construct with the default orientation matrix (identity).
"""
function FEMMHeatDiff(integdomain::ID,
    material::M) where {ID <: IntegDomain, M <: MatHeatDiff}
    return FEMMHeatDiff(integdomain, CSys(manifdim(integdomain.fes)), material)
end

function _buffers1(self::FEMMHeatDiff,
    geom::NodalField{GFT},
    temp::NodalField{FT}) where {GFT, FT}
    # Constants
    fes = self.integdomain.fes
    IntT = eltype(temp.dofnums)
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(temp) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)   # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    Kedim = ndn * nne      # dimension of the element matrix
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of Element coordinates
    elmat = fill(zero(FT), Kedim, Kedim) # buffer
    elvec = fill(zero(FT), Kedim) # buffer
    elvecfix = fill(zero(FT), Kedim) # buffer
    dofnums = fill(zero(IntT), Kedim) # buffer
    loc = fill(zero(FT), 1, sdim) # buffer
    J = fill(zero(FT), sdim, mdim) # buffer
    RmTJ = fill(zero(FT), mdim, mdim) # buffer
    gradN = fill(zero(FT), nne, mdim) # buffer
    kappa_bar = fill(zero(FT), mdim, mdim) # buffer
    kappa_bargradNT = fill(zero(FT), mdim, nne) # buffer
    return ecoords,
    dofnums,
    loc,
    J,
    RmTJ,
    gradN,
    kappa_bar,
    kappa_bargradNT,
    elmat,
    elvec,
    elvecfix
end

"""
    conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{GFT},  temp::NodalField{FT}) where {A<:AbstractSysmatAssembler, GFT, FT}

Compute the conductivity matrix.

# Arguments
- `self` = model machine,
- `assembler` = matrix assembler
- `geom` = geometry field,
- `temp` = temperature field
"""
function conductivity(self::FEMMHeatDiff,
    assembler::A,
    geom::NodalField{GFT},
    temp::NodalField{FT}) where {A <: AbstractSysmatAssembler, GFT, FT}
    mdim = manifdim(finite_elements(self))
    kappa_bar = fill(zero(FT), mdim, mdim) # buffer
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    return bilform_diffusion(self, assembler, geom, temp, DataCache(kappa_bar))
end

function conductivity(self::FEMMHeatDiff,
    geom::NodalField{GFT},
    temp::NodalField{FT}) where {GFT, FT}
    assembler = SysmatAssemblerSparseSymm()
    return conductivity(self, assembler, geom, temp)
end

"""
    energy(self::FEMMHeatDiff, geom::NodalField{GFT},  temp::NodalField{FT}) where {GFT, FT}

Compute the "energy" integral over the interior domain.

The "energy" density is the dot product of the gradient of temperature
and the heat flux.

# Arguments
- `self` = model machine,
- `geom` = geometry field,
- `temp` = temperature field
"""
function energy(self::FEMMHeatDiff,
    geom::NodalField{GFT},
    temp::NodalField{FT}) where {GFT, FT}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    # Prepare assembler and buffers
    ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = _buffers1(self,
        geom,
        temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    gradT = fill(0.0, 1, size(gradN, 2))
    fluxT = reshape(deepcopy(gradT), length(gradT), 1)
    energy = 0.0
    # Now loop over all finite elements in the set
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        gathervalues_asvec!(temp, elvec, fes.conn[i])# retrieve element coordinates
        for j in 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            mulCAtB!(RmTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ)
            mulCAtB!(gradT, reshape(elvec, length(elvec), 1), gradN)
            mulCAB!(fluxT, kappa_bar, reshape(gradT, length(gradT), 1))
            energy += dot(vec(gradT), vec(fluxT)) * (Jac * w[j])
        end # Loop over quadrature points
    end
    return energy
end

"""
    inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{GFT}, u::NodalField{T}, temp::NodalField{FT}, felist::VecOrMat{IntT}, inspector::F, idat, quantity=:heatflux; context...) where {T<:Number, GFT, FT, IntT, F<:Function}

Inspect integration point quantities.

# Arguments
- `geom` - reference geometry field
- `u` - displacement field (ignored)
- `temp` - temperature field
- `felist` - indexes of the finite elements that are to be inspected:
    The fes to be included are: `fes[felist]`.
- `context`    - structure: see the update!() method of the material.
- `inspector` - function with the signature
        `idat = inspector(idat, j, conn, x, out, loc);`
   where
    `idat` - a structure or an array that the inspector may
           use to maintain some state,  for instance heat flux, `j` is the
          element number, `conn` is the element connectivity, `out` is the
          output of the `update!()` method,  `loc` is the location of the
          integration point in the *reference* configuration.
# Output
The updated inspector data is returned.
"""
function inspectintegpoints(self::FEMMHeatDiff,
    geom::NodalField{GFT},
    u::NodalField{T},
    temp::NodalField{FT},
    felist::VecOrMat{IntT},
    inspector::F,
    idat,
    quantity = :heatflux;
    context...) where {T <: Number, GFT, FT, IntT, F <: Function}
    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = _buffers1(self,
        geom,
        temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    # Sort out  the output requirements
    outputcsys = self.mcsys # default: report the vector quantities in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t = 0.0
    dt = 0.0
    Te = fill(zero(FT), nodesperelem(fes)) # nodal temperatures -- buffer
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    qpgradT = fill(zero(FT), 1, sdim) # Temperature gradient -- buffer
    qpflux = fill(zero(FT), sdim) # thermal strain -- buffer
    out1 = fill(zero(FT), sdim) # output -- buffer
    out = fill(zero(FT), sdim)# output -- buffer
    # Loop over  all the elements and all the quadrature points within them
    for ilist in 1:length(felist) # Loop over elements
        i = felist[ilist]
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        gathervalues_asvec!(temp, Te, fes.conn[i])# retrieve element temperatures
        for j in 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            mulCAtB!(RmTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ)
            # Quadrature point quantities
            mulCAB!(qpgradT, reshape(Te, 1, :), gradN) # temperature gradient in material coordinates
            out = update!(self.material,
                qpflux,
                out,
                vec(qpgradT),
                0.0,
                0.0,
                loc,
                fes.label[i],
                quantity)
            if (quantity == :heatflux)   # Transform heat flux vector,  if that is "out"
                mulCAB!(out1, transpose(csmat(self.mcsys)), out)# To global coord sys
                mulCAB!(out, csmat(outputcsys), out1)# To output coord sys
            end
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], ecoords, out, loc)
        end # Loop over quadrature points
    end # Loop over elements
    return idat # return the updated inspector data
end

"""
    capacity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{GFT},  temp::NodalField{FT}) where {A<:AbstractSysmatAssembler, GFT, FT}

Compute the capacity matrix.

# Arguments
- `self` = model machine,
- `assembler` = matrix assembler
- `geom` = geometry field,
- `temp` = temperature field
"""
function capacity(self::FEMMHeatDiff,
    assembler::A,
    geom::NodalField{GFT},
    temp::NodalField{FT}) where {A <: AbstractSysmatAssembler, GFT, FT}
    return bilform_dot(self, assembler, geom, temp, DataCache(self.material.specific_heat))
end

function capacity(self::FEMMHeatDiff,
    geom::NodalField{GFT},
    temp::NodalField{FT}) where {GFT, FT}
    assembler = SysmatAssemblerSparseSymm()
    return capacity(self, assembler, geom, temp)
end

end
