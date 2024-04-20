"""
Heat conduction example described by Amuthan A. Ramabathiran
http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
Unit cube, with known temperature distribution along the boundary,
and uniform heat generation rate inside.
Version: 03/03/2023
"""

module Poisson_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtools.AssemblyModule
using FinEtools.MeshExportModule
using FinEtoolsHeatDiff
using LinearAlgebra
using DataDrop
using SparseArrays
using SymRCM

function Poisson_FE_H20_example(N = 25)
    println("""
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit cube, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular quadratic HEXAHEDRA,
    in a grid of $(N) x $(N) x $(N) edges.
    Version: 02/15/2024
    """)
    t0 = time()
    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3] # conductivity matrix
    Q = -6.0 # internal heat generation rate
    function getsource!(forceout, XYZ, tangents, feid, qpid)
        forceout[1] = Q #heat source
    end
    tempf(x) = (1.0 .+ x[:, 1] .^ 2 + 2.0 .* x[:, 2] .^ 2)#the exact distribution of temperature1
    println("Mesh generation")
    @time fens, fes = H20block(A, A, A, N, N, N)
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))
    println("Searching nodes  for BC")
    Tolerance = 1.0 / N / 100.0
    l1 = selectnode(fens; box = [0.0 0.0 0.0 A 0.0 A], inflate = Tolerance)
    l2 = selectnode(fens; box = [A A 0.0 A 0.0 A], inflate = Tolerance)
    l3 = selectnode(fens; box = [0.0 A 0.0 0.0 0.0 A], inflate = Tolerance)
    l4 = selectnode(fens; box = [0.0 A A A 0.0 A], inflate = Tolerance)
    l5 = selectnode(fens; box = [0.0 A 0.0 A 0.0 0.0], inflate = Tolerance)
    l6 = selectnode(fens; box = [0.0 A 0.0 A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    setebc!(Temp, List, true, 1, tempf(geom.values[List, :])[:])

    numberdofs!(Temp)
    println("Number of free degrees of freedom: $(nfreedofs(Temp))")
    t1 = time()
    material = MatHeatDiff(thermal_conductivity)
    femm = FEMMHeatDiff(IntegDomain(fes, GaussRule(3, 3)), material)
    println("Conductivity")
    @time K = conductivity(femm, geom, Temp)
    println("Internal heat generation")
    # fi = ForceIntensity(Float64, getsource!);# alternative  specification
    fi = ForceIntensity(Float64[Q])
    @time F1 = distribloads(femm, geom, Temp, fi, 3)
    println("Solution of the system")
    @time solve_blocked!(Temp, K, F1)

    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    Error = 0.0
    for k in axes(fens.xyz, 1)
        Error = Error + abs.(Temp.values[k, 1] - tempf(reshape(fens.xyz[k, :], (1, 3)))[1])
    end
    println("Error =$Error")
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    true
end # Poisson_FE_H20_example

function Poisson_FE_T10_example(N = 25)
    println("""
        Mesh of regular QUADRATIC TETRAHEDRA, in a grid of $N x $N x $N edges.
        Version: 03/03/2023
        """)

    t0 = time()

    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3] # conductivity matrix
    Q = -6.0 # internal heat generation rate
    function getsource!(forceout, XYZ, tangents, feid, qpid)
        forceout[1] = Q #heat source
    end
    tempf(x) = (1.0 .+ x[:, 1] .^ 2 + 2.0 .* x[:, 2] .^ 2)#the exact distribution of temperature

    println("Mesh generation")
    fens, fes = T10block(A, A, A, N, N, N)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    println("Searching nodes  for BC")
    Tolerance = 1.0 / N / 100.0
    l1 = selectnode(fens; box = [0.0 0.0 0.0 A 0.0 A], inflate = Tolerance)
    l2 = selectnode(fens; box = [A A 0.0 A 0.0 A], inflate = Tolerance)
    l3 = selectnode(fens; box = [0.0 A 0.0 0.0 0.0 A], inflate = Tolerance)
    l4 = selectnode(fens; box = [0.0 A A A 0.0 A], inflate = Tolerance)
    l5 = selectnode(fens; box = [0.0 A 0.0 A 0.0 0.0], inflate = Tolerance)
    l6 = selectnode(fens; box = [0.0 A 0.0 A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    setebc!(Temp, List, true, 1, tempf(geom.values[List, :])[:])
    applyebc!(Temp)
    numberdofs!(Temp)

    println("Number of free degrees of freedom: $(nfreedofs(Temp))")
    t1 = time()

    material = MatHeatDiff(thermal_conductivity)

    femm = FEMMHeatDiff(IntegDomain(fes, TetRule(4)), material)

    println("Conductivity")
    K = conductivity(femm, geom, Temp)

    println("Internal heat generation")
    # fi = ForceIntensity(Float64, getsource!);# alternative  specification
    fi = ForceIntensity(Float64[Q])
    F1 = distribloads(femm, geom, Temp, fi, 3)
    println("Solution of the system")
    @time solve_blocked!(Temp, K, F1)

    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")

    Error = 0.0
    for k in axes(fens.xyz, 1)
        Error = Error + abs.(Temp.values[k, 1] - tempf(reshape(fens.xyz[k, :], (1, 3)))[1])
    end
    println("Error =$Error")

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

    true
end # Poisson_FE_T10_example

function Poisson_FE_T4_example(N = 25)
    println("""
        Mesh of regular LINEAR TETRAHEDRA, in a grid of $N x $N x $N edges.
        Version: 03/03/2023
        """)
    t0 = time()

    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3] # conductivity matrix
    Q = -6.0 # internal heat generation rate
    function getsource!(forceout, XYZ, tangents, feid, qpid)
        forceout[1] = Q #heat source
    end
    tempf(x) = (1.0 .+ x[:, 1] .^ 2 + 2.0 .* x[:, 2] .^ 2)#the exact distribution of temperature

    println("Mesh generation")
    fens, fes = T4block(A, A, A, N, N, N)

    println("""
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit cube, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular LINEAR TETRAHEDRA,
    in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) nodes).
    Version: 03/03/2023
    """)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    println("Searching nodes  for BC")
    Tolerance = 1.0 / N / 100.0
    l1 = selectnode(fens; box = [0.0 0.0 0.0 A 0.0 A], inflate = Tolerance)
    l2 = selectnode(fens; box = [A A 0.0 A 0.0 A], inflate = Tolerance)
    l3 = selectnode(fens; box = [0.0 A 0.0 0.0 0.0 A], inflate = Tolerance)
    l4 = selectnode(fens; box = [0.0 A A A 0.0 A], inflate = Tolerance)
    l5 = selectnode(fens; box = [0.0 A 0.0 A 0.0 0.0], inflate = Tolerance)
    l6 = selectnode(fens; box = [0.0 A 0.0 A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    setebc!(Temp, List, true, 1, tempf(geom.values[List, :])[:])
    applyebc!(Temp)
    numberdofs!(Temp)

    println("Number of free degrees of freedom: $(nfreedofs(Temp))")
    t1 = time()

    material = MatHeatDiff(thermal_conductivity)

    femm = FEMMHeatDiff(IntegDomain(fes, TetRule(1)), material)

    println("Conductivity")
    K = conductivity(femm, geom, Temp)

    println("Internal heat generation")
    # fi = ForceIntensity(Float64, getsource!);# alternative  specification
    fi = ForceIntensity(Float64[Q])
    F1 = distribloads(femm, geom, Temp, fi, 3)

    println("Solution of the system")
    @time solve_blocked!(Temp, K, F1)

    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")

    Error = 0.0
    for k in axes(fens.xyz, 1)
        Error = Error + abs.(Temp.values[k, 1] - tempf(reshape(fens.xyz[k, :], (1, 3)))[1])
    end
    println("Error =$Error")

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

    true
end # Poisson_FE_T4_example

function zerooutsparse(S)
    S.nzval .= zero(eltype(S.nzval))
    return S
end

function _updcol!(nzval, i, v, st, fi, rowval)
    for r = st:fi
        if i == rowval[r]
            nzval[r] += v
            return
        end
    end
end

function addtosparsepar(S, I, J, V)
    for t in eachindex(J)
        j = J[t]
        _updcol!(S.nzval, I[t], V[t], S.colptr[j], S.colptr[j+1] - 1, S.rowval)
    end
    return S
end

function addtosparse(S, I, J, V)
    for t in eachindex(J)
        j = J[t]
        i = I[t]
        for r = S.colptr[j]:S.colptr[j+1]-1
            if i == S.rowval[r]
                S.nzval[r] += V[t]
                break
            end
        end
    end
    return S
end

function addtosparse3(S, I, J, V)
    @inbounds for t in eachindex(J)
        j = J[t]
        i = I[t]
        for r = S.colptr[j]:S.colptr[j+1]-1
            if i == S.rowval[r]
                S.nzval[r] += V[t]
            end
        end
    end
    return S
end

function addtosparse2(S, I, J, V)
    pj = zero(eltype(J))
    start, stop = zero(eltype(S.colptr)), zero(eltype(S.colptr))
    @inbounds for t in eachindex(J)
        j = J[t]
        if pj != j
            start, stop = S.colptr[j], S.colptr[j+1] - 1
            pj = j
        end
        i = I[t]
        for r = start:stop
            if i == S.rowval[r]
                S.nzval[r] += V[t]
            end
        end
    end
    return S
end

function Poisson_FE_T4_test_sparse_update_example(N = 25)
    println("""
        Mesh of regular LINEAR TETRAHEDRA, in a grid of $N x $N x $N edges.
        Version: 03/03/2023
        """)
    t0 = time()

    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3] # conductivity matrix
    Q = -6.0 # internal heat generation rate
    function getsource!(forceout, XYZ, tangents, feid, qpid)
        forceout[1] = Q #heat source
    end
    tempf(x) = (1.0 .+ x[:, 1] .^ 2 + 2.0 .* x[:, 2] .^ 2)#the exact distribution of temperature

    fens, fes = T4block(A, A, A, N, N, N)
    C = connectionmatrix(FEMMBase(IntegDomain(fes, TetRule(1))), count(fens))
    ordering = symrcm(C)
    fens, fes = reordermesh(fens, fes, ordering)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))

    Tolerance = 1.0 / N / 100.0
    l1 = selectnode(fens; box = [0.0 0.0 0.0 A 0.0 A], inflate = Tolerance)
    l2 = selectnode(fens; box = [A A 0.0 A 0.0 A], inflate = Tolerance)
    l3 = selectnode(fens; box = [0.0 A 0.0 0.0 0.0 A], inflate = Tolerance)
    l4 = selectnode(fens; box = [0.0 A A A 0.0 A], inflate = Tolerance)
    l5 = selectnode(fens; box = [0.0 A 0.0 A 0.0 0.0], inflate = Tolerance)
    l6 = selectnode(fens; box = [0.0 A 0.0 A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    setebc!(Temp, List, true, 1, tempf(geom.values[List, :])[:])
    applyebc!(Temp)
    numberdofs!(Temp)

    material = MatHeatDiff(thermal_conductivity)

    femm = FEMMHeatDiff(IntegDomain(fes, TetRule(1)), material)

    assembler = SysmatAssemblerSparse(0.0)
    setnomatrixresult(assembler, true)
    K = conductivity(femm, assembler, geom, Temp)
    setnomatrixresult(assembler, false)
    @time K = makematrix!(assembler)
    @info "size K = $(size(K)), nnz K = $(nnz(K)), sparsity = $(nnz(K) / prod(size(K)))"
    I, J, V = assembler._rowbuffer, assembler._colbuffer, assembler._matbuffer
    zerooutsparse(K)
    @time addtosparse(K, I, J, V)

    # fi = ForceIntensity(Float64, getsource!);# alternative  specification
    fi = ForceIntensity(Float64[Q])
    F1 = distribloads(femm, geom, Temp, fi, 3)

    println("Solution of the system")
    solve_blocked!(Temp, K, F1)

    Error = 0.0
    for k in axes(fens.xyz, 1)
        Error = Error + abs.(Temp.values[k, 1] - tempf(reshape(fens.xyz[k, :], (1, 3)))[1])
    end
    println("Error =$Error")

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

    true
end # Poisson_FE_T4_example

function _update_buffer_range(elem_mat_nrows, elem_mat_ncols, elem_range, iend)
    buffer_length = elem_mat_nrows * elem_mat_ncols * length(elem_range)
    istart = iend + 1
    iend = iend + buffer_length
    buffer_range = istart:iend
    return buffer_range, iend
end

function _task_local_assembler(a, buffer_range)
    buffer_length = maximum(buffer_range) - minimum(buffer_range) + 1
    matbuffer = view(a._matbuffer, buffer_range)
    rowbuffer = view(a._rowbuffer, buffer_range)
    colbuffer = view(a._colbuffer, buffer_range)
    buffer_pointer = 1
    nomatrixresult = true
    force_init = false
    return SysmatAssemblerSparse(
        buffer_length,
        matbuffer,
        rowbuffer,
        colbuffer,
        buffer_pointer,
        a._row_nalldofs,
        a._col_nalldofs,
        nomatrixresult,
        force_init,
    )
end

function allrun()
    println("#####################################################")
    println("# Poisson_FE_H20_example ")
    Poisson_FE_H20_example()
    println("#####################################################")
    println("# Poisson_FE_T10_example ")
    Poisson_FE_T10_example()
    println("#####################################################")
    println("# Poisson_FE_T4_example ")
    Poisson_FE_T4_example()
    # println("#####################################################")
    # println("# Poisson_FE_T4_altass_example ")
    # Poisson_FE_T4_altass_example()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module Poisson_examples
nothing
