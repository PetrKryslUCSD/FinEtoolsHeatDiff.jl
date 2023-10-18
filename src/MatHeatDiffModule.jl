"""
    MatHeatDiffModule

Module for linear heat diffusion material models.
"""
module MatHeatDiffModule

import FinEtools.MatModule: AbstractMat
using FinEtools.MatrixUtilityModule: mulCAB!

"""
    MatHeatDiff{FT, MTAN<:Function, MUPD<:Function} <: AbstractMat

Type of material model for heat diffusion.
"""
struct MatHeatDiff{FT, MTAN <: Function, MUPD <: Function} <: AbstractMat
    thermal_conductivity::Array{FT, 2}# Thermal conductivity
    specific_heat::FT# Specific heat per unit volume
    mass_density::FT # mass density
    tangentmoduli!::MTAN
    update!::MUPD
end

"""
    MatHeatDiff(thermal_conductivity::Matrix{FT}) where {FT}

Construct material model for heat diffusion.

Supply the matrix of thermal conductivity constants.
"""
function MatHeatDiff(thermal_conductivity::Matrix{FT}) where {FT}
    return MatHeatDiff(thermal_conductivity, zero(FT), zero(FT), tangentmoduli!, update!)
end

"""
    MatHeatDiff(thermal_conductivity::Matrix{FT}, specific_heat::FT) where {FT}

Construct material model for heat diffusion.

Supply the matrix of thermal conductivity constants.
"""
function MatHeatDiff(thermal_conductivity::Matrix{FT}, specific_heat::FT) where {FT}
    return MatHeatDiff(thermal_conductivity,
        specific_heat,
        zero(FT),
        tangentmoduli!,
        update!)
end

"""
    tangentmoduli!(self::MatHeatDiff, kappabar::Matrix{FT}, t = zero(FT), dt = zero(FT), loc::Matrix{FT} = reshape(FT[],0,0), label = 0) where {FT}

Calculate the thermal conductivity matrix.

- `kappabar` = matrix of thermal conductivity (tangent moduli) in material
  coordinate system, supplied as a buffer and overwritten.
"""
function tangentmoduli!(self::MatHeatDiff,
    kappabar::Matrix{FT},
    t = zero(FT),
    dt = zero(FT),
    loc::Matrix{FT} = reshape(FT[], 0, 0),
    label = 0) where {FT}
    copyto!(kappabar, self.thermal_conductivity)
    return kappabar
end

"""
    update!(self::MatHeatDiff, heatflux::Vector{FT}, output::Vector{FT}, gradT::Vector{FT}, t= zero(FT), dt= zero(FT), loc::Matrix{FT}=reshape(FT[],0,0), label=0, quantity=:nothing) where {FT}

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
function update!(self::MatHeatDiff,
    heatflux::Vector{FT},
    output::Vector{FT},
    gradT::Vector{FT},
    t = zero(FT),
    dt = zero(FT),
    loc::Matrix{FT} = reshape(FT[], 0, 0),
    label = 0,
    quantity = :nothing) where {FT}
    sdim = size(self.thermal_conductivity, 2)
    @assert length(heatflux) == sdim
    mulCAB!(heatflux, self.thermal_conductivity, -gradT)
    if quantity == :nothing
        #Nothing to be copied to the output array
    elseif quantity == :heatflux
        (length(output) >= sdim) || (output = zeros(sdim)) # make sure we can store it
        copyto!(output, heatflux)
    end
    return output
end

end
