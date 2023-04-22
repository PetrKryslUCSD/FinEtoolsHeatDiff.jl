"""
    FEMMHeatDiffModule

Module for operations on interiors of domains to construct system
matrices and system vectors for linear heat conduction/diffusion.
"""
module FEMMHeatDiffModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, nodesperelem, manifdim, gradN!
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.CSysModule: CSys, updatecsmat!, csmat
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField
import FinEtools.ElementalFieldModule: ElementalField
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import FinEtools.ForceIntensityModule: ForceIntensity
import FinEtools.FEMMBaseModule: AbstractFEMM, inspectintegpoints
import FinEtools.MatrixUtilityModule: add_gkgt_ut_only!, complete_lt!, locjac!, add_nnt_ut_only!, mulCAtB!, mulCAB!
import LinearAlgebra: norm, dot
import FinEtoolsHeatDiff.MatHeatDiffModule: MatHeatDiff, tangentmoduli!, update!

"""
    FEMMHeatDiff{S<:AbstractFESet, F<:Function, M<:MatHeatDiff} <: AbstractFEMM

Type for heat diffusion finite element modeling machine.
"""
mutable struct FEMMHeatDiff{S<:AbstractFESet, F<:Function, M<:MatHeatDiff} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
end

"""
    FEMMHeatDiff(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M<:MatHeatDiff}

Construct with the default orientation matrix (identity).
"""
function FEMMHeatDiff(integdomain::IntegDomain{S, F}, material::M) where {S<:AbstractFESet, F<:Function, M<:MatHeatDiff}
    return FEMMHeatDiff(integdomain, CSys(manifdim(integdomain.fes)), material)
end

function  _buffers1(self::FEMMHeatDiff, geom::NodalField{FFlt}, temp::NodalField{FFlt})
    # Constants
    fes = self.integdomain.fes
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);   # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    Kedim = ndn*nne;      # dimension of the element matrix
    ecoords = fill(zero(FFlt), nne, ndofs(geom)); # array of Element coordinates
    elmat = fill(zero(FFlt), Kedim, Kedim); # buffer
    elvec = fill(zero(FFlt), Kedim); # buffer
    elvecfix = fill(zero(FFlt), Kedim); # buffer
    dofnums = fill(zero(FInt), Kedim); # buffer
    loc = fill(zero(FFlt), 1, sdim); # buffer
    J = fill(zero(FFlt), sdim, mdim); # buffer
    RmTJ = fill(zero(FFlt), mdim, mdim); # buffer
    gradN = fill(zero(FFlt), nne, mdim); # buffer
    kappa_bar = fill(zero(FFlt), mdim, mdim); # buffer
    kappa_bargradNT = fill(zero(FFlt), mdim, nne); # buffer
    return ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix
end

function  _buffers2(self::FEMMHeatDiff, geom::NodalField{FFlt}, temp::NodalField{FFlt})
    # Constants
    fes = self.integdomain.fes
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);   # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    Kedim = ndn*nne;      # dimension of the element matrix
    ecoords = fill(zero(FFlt), nne, ndofs(geom)); # array of Element coordinates
    elmat = fill(zero(FFlt), Kedim, Kedim); # buffer
    elvec = fill(zero(FFlt), Kedim); # buffer
    elvecfix = fill(zero(FFlt), Kedim); # buffer
    dofnums = fill(zero(FInt), Kedim); # buffer
    loc = fill(zero(FFlt), 1, sdim); # buffer
    J = fill(zero(FFlt), sdim, mdim); # buffer
    # RmTJ = fill(zero(FFlt), mdim, mdim); # buffer
    # gradN = fill(zero(FFlt), nne, mdim); # buffer
    # kappa_bar = fill(zero(FFlt), mdim, mdim); # buffer
    # kappa_bargradNT = fill(zero(FFlt), mdim, nne); # buffer
    return ecoords, dofnums, loc, J, elmat, elvec, elvecfix
end

"""
    conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}

Compute the conductivity matrix.

# Arguments
- `self` = model machine,
- `assembler` = matrix assembler
- `geom` = geometry field,
- `temp` = temperature field
"""
function conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and buffers
    ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat = _buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    startassembly!(assembler, size(elmat)..., count(fes), temp.nfreedofs, temp.nfreedofs);
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        fill!(elmat,  0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            mulCAtB!(RmTJ,  csmat(self.mcsys),  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ);
            # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
            add_gkgt_ut_only!(elmat, gradN, (Jac*w[j]), kappa_bar, kappa_bargradNT)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(temp, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function conductivity(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt})
    assembler = SysmatAssemblerSparseSymm();
    return conductivity(self, assembler, geom, temp);
end

"""
    nzebcloadsconductivity(self::FEMMHeatDiff, assembler::A,  geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}

Compute load vector for nonzero EBC of prescribed temperature.

# Arguments
- `self` = model machine,
- `assembler` = matrix assembler
- `geom` = geometry field,
- `temp` = temperature field
"""
function nzebcloadsconductivity(self::FEMMHeatDiff, assembler::A,  geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and buffers
    ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = _buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    startassembly!(assembler, temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        gatherfixedvalues_asvec!(temp, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) != 0. # Is the load nonzero?
            fill!(elmat,  0.0);
            for j=1:npts # Loop over quadrature points
                locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                updatecsmat!(self.mcsys, loc, J, fes.label[i]);
                mulCAtB!(RmTJ, csmat(self.mcsys), J); # local Jacobian matrix
                gradN!(fes, gradN, gradNparams[j], RmTJ);
                # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
                add_gkgt_ut_only!(elmat, gradN, (Jac*w[j]), kappa_bar, kappa_bargradNT)
            end # Loop over quadrature points
            complete_lt!(elmat)
            mulCAB!(elvec, elmat, elvecfix) # compute  the load vector
            gatherdofnums!(temp, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
        end
    end
    return makevector!(assembler);
end

function nzebcloadsconductivity(self::FEMMHeatDiff,  geom::NodalField{FFlt},   temp::NodalField{FFlt})
    assembler = SysvecAssembler()
    return  nzebcloadsconductivity(self, assembler, geom, temp);
end

"""
    energy(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt})

Compute the "energy" integral over the interior domain.

The "energy" density is the dot product of the gradient of temperature
and the heat flux.

# Arguments
- `self` = model machine,
- `geom` = geometry field,
- `temp` = temperature field
"""
function energy(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt})
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and buffers
    ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = _buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    gradT = fill(0.0, 1, size(gradN, 2))
    fluxT = reshape(deepcopy(gradT), length(gradT), 1)
    energy = 0.0
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        gathervalues_asvec!(temp, elvec, fes.conn[i]);# retrieve element coordinates
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            mulCAtB!(RmTJ,  csmat(self.mcsys),  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ);
            mulCAtB!(gradT, reshape(elvec, length(elvec), 1), gradN)
           	mulCAB!(fluxT, kappa_bar, reshape(gradT, length(gradT), 1))
            energy += dot(vec(gradT), vec(fluxT)) * (Jac*w[j])
        end # Loop over quadrature points
    end
    return energy;
end

"""
    inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{FFlt}, u::NodalField{T}, temp::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:heatflux; context...) where {T<:Number, F<:Function}

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
function inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{FFlt}, u::NodalField{T}, temp::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:heatflux; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    ecoords, dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = _buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = tangentmoduli!(self.material, kappa_bar)
    # Sort out  the output requirements
    outputcsys = self.mcsys; # default: report the vector quantities in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t= 0.0
    dt = 0.0
    Te = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    qpgradT = fill(zero(FFlt), 1, sdim); # Temperature gradient -- buffer
    qpflux = fill(zero(FFlt), sdim); # thermal strain -- buffer
    out1 = fill(zero(FFlt), sdim); # output -- buffer
    out =  fill(zero(FFlt), sdim);# output -- buffer
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist];
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        gathervalues_asvec!(temp, Te, fes.conn[i]);# retrieve element temperatures
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            mulCAtB!(RmTJ,  csmat(self.mcsys),  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ);
            # Quadrature point quantities
            mulCAB!(qpgradT, reshape(Te, 1, :), gradN); # temperature gradient in material coordinates
            out = update!(self.material, qpflux, out, vec(qpgradT), 0.0, 0.0, loc, fes.label[i], quantity)
            if (quantity == :heatflux)   # Transform heat flux vector,  if that is "out"
                mulCAB!(out1, transpose(csmat(self.mcsys)), out);# To global coord sys
                mulCAB!(out, csmat(outputcsys), out1);# To output coord sys
            end
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], ecoords, out, loc);
        end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

"""
    capacity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}

Compute the capacity matrix.

# Arguments
- `self` = model machine,
- `assembler` = matrix assembler
- `geom` = geometry field,
- `temp` = temperature field
"""
function capacity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}
	fes = self.integdomain.fes
	ecoords, dofnums, loc, J, elmat, elvec, elvecfix = _buffers2(self, geom, temp)
	# Precompute basis f. values + basis f. gradients wrt parametric coor
	npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
	# Material
	specific_heat  =  self.material.specific_heat;
	startassembly!(assembler, size(elmat)..., count(fes), temp.nfreedofs, temp.nfreedofs);
	for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
		fill!(elmat, 0.0); # Initialize element matrix
		for j = 1:npts # Loop over quadrature points
			locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
			Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
			ffactor = Jac*specific_heat*w[j]
			add_nnt_ut_only!(elmat, Ns[j], ffactor)
		end # Loop over quadrature points
		complete_lt!(elmat)
		gatherdofnums!(temp, dofnums, fes.conn[i]);# retrieve degrees of freedom
		assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
	end # Loop over elements
	return makematrix!(assembler);
end

function capacity(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt})
    assembler = SysmatAssemblerSparseSymm();
    return capacity(self, assembler, geom, temp);
end

end
