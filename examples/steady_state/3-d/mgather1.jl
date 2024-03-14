module mgather1
using Test
using Random

function test()
    N = 10000 # Number of nodes in the mesh
    nen = 20 # Number of nodes per element
    nloops = 2 * N
    all_indexes = [randperm(N)[1:nen] for _ in 1:nloops]
    buffnen3 = rand(nen, 3)
    buff3nen = rand(3, nen)
    dataN3 = rand(N, 3)
    data3N = Matrix(transpose(dataN3))

    t1 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for i in 1:nen
            ii = indexes[i]
            for j in 1:3
                buffnen3[i, j] = dataN3[ii, j]
            end
        end
    end

    t2 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for j in 1:3
          # alternative to the loop below 
          # buffnen3[:, j] .= dataN3[indexes, j] # SLOW
            @inbounds for i in 1:nen
                ii = indexes[i]
                buffnen3[i, j] = dataN3[ii, j]
            end
        end
    end

    t3 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for i in 1:nen
            ii = indexes[i]
            for j in 1:3
                buff3nen[j, i] = dataN3[ii, j]
            end
        end
    end

    t4 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for j in 1:3
            for i in 1:nen
                ii = indexes[i]
                buff3nen[j, i] = dataN3[ii, j]
            end
        end
    end


    t5 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for i in 1:nen
            ii = indexes[i]
            for j in 1:3
                buffnen3[i, j] = data3N[j, ii]
            end
        end
    end

    t6 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for j in 1:3
            for i in 1:nen
                ii = indexes[i]
                buffnen3[i, j] = data3N[j, ii]
            end
        end
    end

    t7 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for i in 1:nen
            ii = indexes[i]
            for j in 1:3
                buff3nen[j, i] = data3N[j, ii]
            end
        end
    end

    t8 = @elapsed for l = 1:nloops
        indexes = view(all_indexes, l)[1]
        @inbounds for j in 1:3
            for i in 1:nen
                ii = indexes[i]
                buff3nen[j, i] = data3N[j, ii]
            end
        end
    end

    [t1, t2, t3, t4, t5, t6, t7, t8] ./ nloops .* 1e6 # In microseconds
end
end
using Main.mgather1
ts = [0.0 for i in 1:8]
ntries = 10
for i in 1:ntries
  @info "Try $i"
    ts .+= mgather1.test()
end
ts ./= ntries
ts = Float32.(ts)

println("Mesh data N x 3, Element buffer nen x 3, Loop i, j: Time $(ts[1]) [mus]")
println("Mesh data N x 3, Element buffer nen x 3, Loop j, i: Time $(ts[2]) [mus]")
println("Mesh data N x 3, Element buffer 3 x nen, Loop i, j: Time $(ts[3]) [mus]")
println("Mesh data N x 3, Element buffer 3 x nen, Loop j, i: Time $(ts[4]) [mus]")
println("Mesh data 3 x N, Element buffer nen x 3, Loop i, j: Time $(ts[5]) [mus]")
println("Mesh data 3 x N, Element buffer nen x 3, Loop j, i: Time $(ts[6]) [mus]")
println("Mesh data 3 x N, Element buffer 3 x nen, Loop i, j: Time $(ts[7]) [mus]")
println("Mesh data 3 x N, Element buffer 3 x nen, Loop j, i: Time $(ts[8]) [mus]")
