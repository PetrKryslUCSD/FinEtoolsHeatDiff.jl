"""
Transient heating of a rectangular steel block suddenly subject to surface temperature of 100° C.

The initial temperature of the block is 20° C. The temperature at the center
climbs from 20° C to approximately 90° C in 1200 seconds.
"""

module steel_block_examples
using FinEtools
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
import LinearAlgebra: cholesky, mul!
using PlotlyLight

function steel_block_m()
    # Thermal conductivity
    thermal_conductivity = [i == j ? 44.0 * phun("W/K/m") : zero(FFlt)
                            for i in 1:3, j in 1:3] # conductivity matrix
    # Specific heat: the specific heat is per unit volume
    rho = 8000 * phun("kg/m^3")
    specific_heat = 470.0 * phun("J/kg/K") * rho
    theta = 1.0 # generalized trapezoidal method parameter
    dt = 1.0 # time step
    tend = 1200 # length of the time interval
    T0(x, y, z) = 20.0# Initial distribution of temperature
    Ts = 100.0 # Prescribed surface temperature
    dimensions = (0.2, 0.15, 0.3) .* phun("m")
    nels = (10, 10, 10)

    tolerance = min(dimensions...) / max(nels...) / 100

    fens, fes = H8block(dimensions..., nels...)
    @show lcenter = selectnode(fens; nearestto = [0.0 0.0 0.0], inflate = tolerance)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    # Initial temperature
    for i in 1:count(fens)
        Temp.values[i] = T0(fens.xyz[i, :]...)
    end

    lx = selectnode(fens;
        box = [dimensions[1] dimensions[1] -Inf Inf -Inf Inf],
        inflate = tolerance)
    ly = selectnode(fens;
        box = [-Inf Inf dimensions[2] dimensions[2] -Inf Inf],
        inflate = tolerance)
    lz = selectnode(fens;
        box = [-Inf Inf -Inf Inf dimensions[3] dimensions[3]],
        inflate = tolerance)
    for (i, l) in enumerate([lx, ly, lz])
        setebc!(Temp, l, true, 1, Ts)
    end
    applyebc!(Temp)
    numberdofs!(Temp)

    centerdof = Temp.dofnums[lcenter]

    material = MatHeatDiff(thermal_conductivity, specific_heat)
    femm = FEMMHeatDiff(IntegDomain(fes, GaussRule(3, 3)), material)

    K = conductivity(femm, geom, Temp)
    C = capacity(femm, geom, Temp)
    L = nzebcloadsconductivity(femm, geom, Temp)

    A = cholesky((1 / dt * C + theta * K))
    B = (1 / dt * C - (1 - theta) * K)

    Tn = gathersysvec(Temp)
    Tn1 = deepcopy(Tn)

    tstart = time()
    t = 0.0
    ts = Float64[t]
    @show Tn1[centerdof, 1][1]
    center_T = Float64[Tn1[centerdof, 1][1]]
    while true # Time stepping
        t = t + dt # Advance the current time
        Tn .= Tn1 # current temperature system vector

        # Solve for the new temperature  vector
        rhs = B * Tn + L # complete right-hand side (loads)
        Tn1 .= A \ rhs # compute next temperature
        push!(ts, t)
        push!(center_T, Tn1[centerdof, 1][1])
        if (t == tend) # Have we reached the end?  If so jump out.
            break
        end
        if (t + dt > tend) # Adjust the last time step so that we exactly reach tend
            dt = tend - t
        end
    end
    println("Time = $(time() - tstart)")

    @show minimum(center_T), maximum(center_T)

    plt = PlotlyLight.Plot()
    plt(x = vec(ts), y = vec(center_T))
    plt.layout.xaxis.title = "t"
    plt.layout.yaxis.title = "T"
    display(plt)
    # @gp "set terminal windows 0 " :-
    # @gp  :- vec(geom.values) vec(Temp.values) "lw 2 lc rgb 'red' with lines title 'Temperature at the corner' "
    # @gp  :- "set xlabel 'Time [s]'"
    # @gp  :- "set ylabel 'Temperature [degrees C]'"

    true
end

function allrun()
    println("#####################################################")
    println("# steel_block_m ")
    steel_block_m()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module steel_block
nothing
