module rectangle_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
using Targe2
using Targe2.Utilities: polygon

function rectangle_1_T3_example()
    println("""
    rectangle 1
    """)
    kappa = [1.0 0; 0 1.0] # conductivity matrix
    a, d = 1.0, 1.6
    v = [0 0; a 0; a d; 0 d]

    thickness = 1.0
    tolerance = min(a, d) / 10000

    s = polygon(v)

    commands = """
   $([_s * "\n" for _s in s]...)
   subregion 1  property 1 boundary $([string(_i) * " " for _i in eachindex(s)]...)
   m-ctl-point constant 0.03
   """

    mesh = triangulate(commands)
    fens = FENodeSet(mesh.xy)
    fes = FESetT3(mesh.triconn)

    edge_fes = meshboundary(fes)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    l1 = connectednodes(edge_fes)
    setebc!(Temp, l1, 1, zero(Float64))
    applyebc!(Temp)

    numberdofs!(Temp)

    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1)), material)

    K = conductivity(femm, geom, Temp)

    fi = ForceIntensity(Float64[1.0])
    F = distribloads(femm, geom, Temp, fi, 3)

    solve_blocked!(Temp, K, F)

    File = "rectangle_1_T3_example-T.vtk"
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

    File = "rectangle_1_T3_example-q.vtk"
    vtkexportvectors(File, qplocs, [("heatflux", qpfluxes)])

    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

    true
end # rectangle_1_T3_example

function allrun()
    println("#####################################################")
    println("# rectangle_1_T3_example ")
    rectangle_1_T3_example()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module rectangle_examples
nothing
