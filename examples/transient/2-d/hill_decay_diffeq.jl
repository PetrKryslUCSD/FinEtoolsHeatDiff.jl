module hill_decay_diffeq
using FinEtools
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule    
import LinearAlgebra: cholesky
using UnicodePlots
using DifferentialEquations

function hill_decay_t3()
	thermal_conductivity =  [i==j ? 0.2 : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
	Width = 60.0
	Height = 40.0
	N = 50
	specific_heat = 1.0
	theta = 1/2; # generalized trapezoidal method parameter
	dt = 2.0; # time step
	tend = 200*dt; # length of the time interval
	T0(x, y) = 500.0/(x^2+y^2+5);# Initial distribution of temperature
	
	tolerance =Width/N/100;

    fens,fes =T3block(Width, Height, N, N)
    for i in 1:count(fens)
    	fens.xyz[i, 1] -= Width / 2
    	fens.xyz[i, 2] -= Height / 2
    end
    l1  = selectnode(fens; box=[Width/2 Width/2 Height/2 Height/2],  inflate = tolerance)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))

    for i in 1:count(fens)
    	Temp.values[i] = T0(fens.xyz[i, :]...)
    end

    applyebc!(Temp)
    numberdofs!(Temp)

    cornerdof = Temp.dofnums[l1]

    material = MatHeatDiff(thermal_conductivity, specific_heat)
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1)), material)

    K = conductivity(femm, geom, Temp)
    C = capacity(femm, geom, Temp)
    
    Tn = gathersysvec(Temp)
    
	f(dT, T, p, t) = begin dT .= -(C \ (K * T)); end
	tspan = (0.0, tend)
	prob = ODEProblem(f, Tn, tspan)
	sol = solve(prob, Rosenbrock23(), dtmin = dt/100)
	@show sol

	# ts = Float64[]
	# Corner_T = Float64[];
	# for t in 0.0:dt:tend 
	# 	T = sol(t)
	# 	push!(ts, t)
	# 	push!(Corner_T, T[cornerdof])
	# end
	
	# @show minimum(Corner_T), maximum(Corner_T)
	# plt = lineplot(vec(ts), vec(Corner_T), canvas = DotCanvas, title = "Transient temperature at the corner", name = "T", xlabel = "Time", ylabel = "T")
	# display(plt)
	# File =  "a.vtk"
	# MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

    true
end

function allrun()
    println("#####################################################")
    println("# hill_decay_t3 ")
    hill_decay_t3()
end # function allrun

end # module hill_decay_examples

using .hill_decay_diffeq
hill_decay_diffeq.hill_decay_t3()