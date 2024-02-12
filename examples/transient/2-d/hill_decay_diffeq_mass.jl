module hill_decay_diffeq_mass
using FinEtools
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
import LinearAlgebra: lu!, mul!, cholesky, I
using UnicodePlots
using DifferentialEquations

function hill_decay_t3()
    thermal_conductivity = [i == j ? 0.2 : zero(FFlt) for i = 1:2, j = 1:2] # conductivity matrix
    Width = 60.0
    Height = 40.0
    N = 50
    specific_heat = 1.0
    dt = 2.0 # time step
    tend = 20 * dt # length of the time interval
    T0(x, y) = 500.0 / (x^2 + y^2 + 5)# Initial distribution of temperature

    tolerance = Width / N / 100

    fens, fes = T3block(Width, Height, N, N)
    for i = 1:count(fens)
        fens.xyz[i, 1] -= Width / 2
        fens.xyz[i, 2] -= Height / 2
    end
    l1 = selectnode(
        fens;
        box = [Width / 2 Width / 2 Height / 2 Height / 2],
        inflate = tolerance,
    )

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    for i = 1:count(fens)
        Temp.values[i] = T0(fens.xyz[i, :]...)
    end

    applyebc!(Temp)
    numberdofs!(Temp)

    cornerdof = Temp.dofnums[l1]

    material = MatHeatDiff(thermal_conductivity, specific_heat)
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(3)), material)

    # The ODE to solve is C*dT/dt + K*T = 0
    K = conductivity(femm, geom, Temp)
    C = capacity(femm, geom, Temp)
    CF = cholesky(C)

    Tn = gathersysvec(Temp)

    tstart = time()
    f(dT, T, p, t) = begin
        dT .= -((K * T))
    end
    tspan = (0.0, tend)
    ff = ODEFunction(f, mass_matrix = C)
    prob = ODEProblem(ff, Tn, tspan)
    sol = solve(prob, Rosenbrock23(autodiff = false), abstol = 1.0e-4, reltol = 1.0e-4)
    # @show sol
    println("Time = $(time() - tstart)")

    # This is postprocessing  to extract the data for the plot.
    ts = Float64[]
    Corner_T = Float64[]
    for t in sol.t
        T = sol(t)
        push!(ts, t)
        push!(Corner_T, T[cornerdof][1])
    end

    @show ts
    @show Corner_T
    @show minimum(Corner_T), maximum(Corner_T)
    plt = lineplot(
        vec(ts),
        vec(Corner_T),
        canvas = DotCanvas,
        title = "Transient temperature at the corner",
        name = "T",
        xlabel = "Time",
        ylabel = "T",
    )
    display(plt)

    true
end

end # module hill_decay_diffeq_mass

using .hill_decay_diffeq_mass
hill_decay_diffeq_mass.hill_decay_t3()
