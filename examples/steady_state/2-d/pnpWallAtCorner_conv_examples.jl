module pnpWallAtCorner_conv_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtools.MeshExportModule
using LinearAlgebra: cholesky

function pnpWallAtCorner_conv()
    println("""
    pnpWallAtCorner_conv.
    Two-dimensional heat transfer with convection.
    """)

    X = Float64[-500 250; 250 250; -500 0; 0 0; -500 300; 300 300]
    ki = 0.05
    kc = 1.8
    Dz = 1.0
    hi = 5.0 / 1000
    he = 20.0 / 1000

    

    fens = FENodeSet(X)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    numberdofs!(Temp)
    Temp.dofnums[:] = 1:size(X, 1)
    ambtemp = deepcopy(Temp)
    ambtemp.values[:] = [0, 0, 20, 20, -10, -10]

    # Insulation Layer
    fes = FESetT3([1 2 5; 2 6 5; ])
    m = MatHeatDiff([ki 0; 0 ki])
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1), Dz), m)
    K = conductivity(femm, geom, Temp)

    # Concrete layer
    fes = FESetT3([1 4 2; 3 4 1; ])
    m = MatHeatDiff([kc 0; 0 kc])
    femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1), Dz), m)
    K += conductivity(femm, geom, Temp)

    bfes = FESetL2([3 4;])
    cfemm = FEMMHeatDiffSurf(IntegDomain(bfes, GaussRule(1, 2), Dz), hi)
    cfemm = FEMMHeatDiffSurf(IntegDomain(bfes, TrapezoidalRule(1), Dz), hi)
    H = surfacetransfer(cfemm, geom, Temp)
    L = surfacetransferloads(cfemm, geom, Temp, ambtemp)

    bfes = FESetL2([6 5;])
    cfemm = FEMMHeatDiffSurf(IntegDomain(bfes, GaussRule(1, 2), Dz), he)
    cfemm = FEMMHeatDiffSurf(IntegDomain(bfes, TrapezoidalRule(1), Dz), he)
    H += surfacetransfer(cfemm, geom, Temp)
    L += surfacetransferloads(cfemm, geom, Temp, ambtemp)
    
    # F3 = surfacetransferloads(femm, geom, temp, amb);
    K_ff, K_fd = matrix_blocked(K, nfreedofs(Temp), nfreedofs(Temp))[(:ff, :fd)]
    @show H_ff, H_fd = matrix_blocked(H, nfreedofs(Temp), nfreedofs(Temp))[(:ff, :fd)]
    L_f = L[:]
    T_f = (K_ff + H_ff) \ (L_f)
    scattersysvec!(Temp, T_f)

    @show Temp
end

function allrun()
    #println("#####################################################")
    println("# pnpWallAtCorner_conv ")
    pnpWallAtCorner_conv()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module pnpConcreteColumnsymb_conv_examples
nothing
