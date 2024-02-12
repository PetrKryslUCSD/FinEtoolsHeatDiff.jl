module maze_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
using Targe2
using Targe2.Utilities: polygon

function maze_1_T3_example()
    println("""
    Maze 1
    """)
    kappa = [1.0 0; 0 1.0] # conductivity matrix
    magn = 0.09# heat flux along the boundary
    a, b, d = 2.3, 0.05, 12
    v = [0 0; a 0; a d; a+b d; a+b 0; 2*a+b 0; 2*a+b d+a; 0 d+a]

    thickness = 1.0
    tolerance = min(a, b, d) / 10000

    s = polygon(v)

    commands = """
   $([_s * "\n" for _s in s]...)
   subregion 1  property 1 boundary $([string(_i) * " " for _i in eachindex(s)]...)
   m-ctl-point constant 0.25
   """

    mesh = triangulate(commands)
    fens = FENodeSet(mesh.xy)
    fes = FESetT3(mesh.triconn)

    edge_fes = meshboundary(fes)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    l1 = selectnode(fens; box = [0.0 0.0 0.0 0.0], inflate = tolerance)
    setebc!(Temp, l1, 1, zero(FFlt))
    applyebc!(Temp)

    numberdofs!(Temp)

    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1)), material)

    K = conductivity(femm, geom, Temp)

    l1 = selectelem(fens, edge_fes, box = [0 a 0.0 0.0])
    el1femm = FEMMBase(IntegDomain(subset(edge_fes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[-magn])#entering the domain
    F1 = (-1) * distribloads(el1femm, geom, Temp, fi, 2)

    l2 = selectelem(fens, edge_fes, box = [a + b a + b + a 0.0 0.0])
    el2femm = FEMMBase(IntegDomain(subset(edge_fes, l2), GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[+magn])#leaving the domain
    F2 = (-1) * distribloads(el2femm, geom, Temp, fi, 2)

    solve_blocked!(Temp, K, F1 + F2)

    File = "maze_1_T3_example-T.vtk"
    vtkexportmesh(
        File,
        connasarray(fes),
        [geom.values Temp.values],
        FinEtools.MeshExportModule.VTK.T3;
        scalars = [("Temperature", Temp.values)],
    )
    @async run(`"paraview.exe" $File`)

    qplocs = []
    qpfluxes = []
    function inspector(idat, i, conn, xe, out, loc)
        qplocs, qpfluxes = idat
        push!(qplocs, copy(loc))
        push!(qpfluxes, copy(out))
        return (qplocs, qpfluxes)
    end
    idat = inspectintegpoints(
        femm,
        geom,
        NodalField([1.0]),
        Temp,
        collect(1:count(fes)),
        inspector,
        (qplocs, qpfluxes),
        :heatflux,
    )

    File = "maze_1_T3_example-q.vtk"
    vtkexportvectors(File, qplocs, [("heatflux", qpfluxes)])

    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

    true
end # maze_1_T3_example

function allrun()
    println("#####################################################")
    println("# maze_1_T3_example ")
    maze_1_T3_example()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module maze_examples
nothing
