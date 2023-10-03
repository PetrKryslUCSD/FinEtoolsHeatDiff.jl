"""
    MatHeatDiffModule

Module for linear heat diffusion material models.
"""
module MatHeatDiffModule

import FinEtools.MatModule: AbstractMat
using FinEtools.MatrixUtilityModule: mulCAB!

"""
    MatHeatDiff{FloatT, MTAN<:Function, MUPD<:Function} <: AbstractMat

Type of material model for heat diffusion.
"""
struct MatHeatDiff{FloatT, MTAN<:Function, MUPD<:Function} <: AbstractMat
	thermal_conductivity::Array{FloatT, 2};# Thermal conductivity
	specific_heat::FloatT;# Specific heat per unit volume
	mass_density::FloatT # mass density
	tangentmoduli!::MTAN
	update!::MUPD
end

"""
    MatHeatDiff(thermal_conductivity::Matrix{FloatT}) where {FloatT}

Construct material model for heat diffusion.

Supply the matrix of thermal conductivity constants.
"""
function MatHeatDiff(thermal_conductivity::Matrix{FloatT}) where {FloatT}
    return MatHeatDiff(thermal_conductivity, zero(FloatT), zero(FloatT), tangentmoduli!, update!)
end

"""
    MatHeatDiff(thermal_conductivity::Matrix{FloatT}, specific_heat::FloatT) where {FloatT}

Construct material model for heat diffusion.

Supply the matrix of thermal conductivity constants.
"""
function MatHeatDiff(thermal_conductivity::Matrix{FloatT}, specific_heat::FloatT) where {FloatT}
    return MatHeatDiff(thermal_conductivity, specific_heat, zero(FloatT), tangentmoduli!, update!)
end

"""
    tangentmoduli!(self::MatHeatDiff, kappabar::Matrix{FloatT}, t = zero(FloatT), dt = zero(FloatT), loc::Matrix{FloatT} = reshape(FloatT[],0,0), label = 0) where {FloatT}

Calculate the thermal conductivity matrix.

- `kappabar` = matrix of thermal conductivity (tangent moduli) in material
  coordinate system, supplied as a buffer and overwritten.
"""
function tangentmoduli!(self::MatHeatDiff, kappabar::Matrix{FloatT}, t = zero(FloatT), dt = zero(FloatT), loc::Matrix{FloatT} = reshape(FloatT[],0,0), label = 0) where {FloatT}
    copyto!(kappabar, self.thermal_conductivity);
    return kappabar
end

"""
    update!(self::MatHeatDiff, heatflux::Vector{FloatT}, output::Vector{FloatT}, gradT::Vector{FloatT}, t= zero(FloatT), dt= zero(FloatT), loc::Matrix{FloatT}=reshape(FloatT[],0,0), label=0, quantity=:nothing) where {FloatT}

Update material state.

# Arguments
- `gradT` = thermal gradient vector,
- `t` = current time,
- `dt` = current time step,
- `loc` = location of the quadrature point in global Cartesian coordinates,
- `label` = label of the finite element in which the quadrature point is located.
- `quantity` = quantity to be output (`:heatflux`)

# Output
- `heatflux` = heat flux vector, allocated by the caller with a size of
  the embedding space. The components of the heat flux vector are
  calculated and stored in the `heatflux` vector.
- `output` =  array which is (if necessary) allocated  in an appropriate size, filled with the output quantity, and returned.
"""
function update!(self::MatHeatDiff, heatflux::Vector{FloatT}, output::Vector{FloatT}, gradT::Vector{FloatT}, t= zero(FloatT), dt= zero(FloatT), loc::Matrix{FloatT}=reshape(FloatT[],0,0), label=0, quantity=:nothing) where {FloatT}
	sdim = size(self.thermal_conductivity, 2)
    @assert length(heatflux) == sdim
    mulCAB!(heatflux, self.thermal_conductivity, -gradT);
    if quantity == :nothing
        #Nothing to be copied to the output array
    elseif quantity == :heatflux
        (length(output) >= sdim) || (output = zeros(sdim)) # make sure we can store it
        copyto!(output, heatflux);
    end
    return output
end

end
