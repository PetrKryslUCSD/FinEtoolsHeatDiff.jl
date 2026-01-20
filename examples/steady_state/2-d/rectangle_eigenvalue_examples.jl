module rectangle_eigenvalue_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using LinearAlgebra
using Targe2
using Targe2.Utilities: polygon
using Arpack
using SparseArrays

function rectangle_1_T3_example(h)
    
    kappa = [1.0 0; 0 1.0] # conductivity matrix
    a, d = 1.0, 1.6
    v = [0 0; a 0; a d; 0 d]

    s = polygon(v)

    commands = """
   $([_s * "\n" for _s in s]...)
   subregion 1  property 1 boundary $([string(_i) * " " for _i in eachindex(s)]...)
   m-ctl-point constant $(h)
   """

    mesh = triangulate(commands)
    fens = FENodeSet(mesh.xy)
    fes = FESetT3(mesh.triconn)
    println("""
    Rectangle h=$(h), #el=$(count(fes))
    """)

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

    D = spdiagm(1.0 ./ sqrt.(diag(K)))
    K = D * K * D
    
    fi = ForceIntensity(Float64[1.0])
    F = distribloads(femm, geom, Temp, fi, 3)
    F = D * F

    solve_blocked!(Temp, K, F)
    Tbar = gathersysvec(Temp, DOF_KIND_ALL)
    T = D  \ Tbar
    scattersysvec!(Temp, T, DOF_KIND_ALL)

    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
    
    Kff = matrix_blocked_ff(K, nfreedofs(Temp))
    eval, evec = eigs(Symmetric(h^-2 * Kff); nev=3, which=:SM)
    @show eval 

    File = "rectangle_1_T3_example-h=$(round(h, digits=5))-T.vtk"
    scalars = Any[("Temperature", deepcopy(Temp.values))]
    for (i, e) in enumerate(eval)
        scattersysvec!(Temp, evec[:, i])
        push!(scalars, ("ev$i", deepcopy(Temp.values)))
    end
    vtkexportmesh(
        File,
        connasarray(fes),
        [geom.values Temp.values],
        FinEtools.MeshExportModule.VTK.T3;
        scalars = scalars,
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

    File = "rectangle_1_T3_example-h=$(round(h, digits=5))-q.vtk"
    vtkexportvectors(File, qplocs, [("heatflux", qpfluxes)])

    true
end # rectangle_1_T3_example

function eigenvalues()
    # for h in 0.03511 ./ (1:2:3)
    for h in 0.0511 ./ [2^i for i in 0:1:4]
        rectangle_1_T3_example(h)
    end
end

function allrun()
    println("#####################################################")
    println("# eigenvalues ")
    eigenvalues()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module rectangle_examples
nothing
