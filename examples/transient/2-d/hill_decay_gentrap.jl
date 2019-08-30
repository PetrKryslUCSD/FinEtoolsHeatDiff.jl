module hill_decay_gentrap
using FinEtools
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule    
import LinearAlgebra: cholesky
using UnicodePlots

function hill_decay_t3()
	thermal_conductivity =  [i==j ? 0.2 : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
	Width = 60.0
	Height = 40.0
	N = 10
	specific_heat = 1.0
	theta = 1.0; # generalized trapezoidal method parameter
	dt = 2.0; # time step
	tend = 20*dt; # length of the time interval
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
    
    A = cholesky((1/dt*C+theta*K))
    B = (1/dt*C-(1-theta)*K)
    
    Tn = gathersysvec(Temp)
    Tn1 = deepcopy(Tn)

    t = 0.0
    ts = Float64[t]
    @show Tn1[cornerdof, 1][1]
    Corner_T = Float64[Tn1[cornerdof, 1][1]];
    while true # Time stepping
    	t = t+dt; # Advance the current time
    	Tn .= Tn1 # current temperature system vector

    	# Solve for the new temperature  vector
    	rhs = B*Tn; # complete right-hand side (loads)
    	Tn1 .= A\rhs; # compute next temperature
    	push!(ts, t)
    	push!(Corner_T, Tn1[cornerdof, 1][1])
    	if (t == tend) # Have we reached the end?  If so jump out.
	    	break;
	    end
	    if (t+dt > tend) # Adjust the last time step so that we exactly reach tend
	    	dt = tend-t;
	    end
	end
     
     @show minimum(Corner_T), maximum(Corner_T)
     plt = lineplot(vec(ts), vec(Corner_T), canvas = DotCanvas, title = "Transient temperature at the corner", name = "T", xlabel = "Time", ylabel = "T")
     display(plt)
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

    true
end

end # module hill_decay_gentrap

using .hill_decay_gentrap
hill_decay_gentrap.hill_decay_t3()
