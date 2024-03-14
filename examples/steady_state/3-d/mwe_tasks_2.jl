module mwe_tasks_2
using Base.Threads
function work(r)
    s = 0.0
    for j in r
        s = s + exp(-(j - minimum(r))^2 / (maximum(r) - minimum(r))^2)
    end
    s
end
function test(nchunks = 2)
    # @info "nchunks = $(nchunks)"
    # @info "nthreads.((:default, :interactive)) = $(Threads.nthreads.((:default, :interactive)))"
    # @info "maxthreadid = $(Threads.maxthreadid())"
    # @info "threadpool.(1:Threads.maxthreadid()) = $(threadpool.(1:Threads.maxthreadid()))"
    # @info "Threads.threadpoolsize.((:default, :interactive)) = $(Threads.threadpoolsize.((:default, :interactive)))"
    N = 200000000
    chincr = N / nchunks
    @assert nchunks * chincr == N
    chunks = [(((i-1)*chincr+1:i*chincr), i) for i = 1:nchunks]

    s = Float64[]
    start = time()
    Threads.@sync begin
        for ch in chunks
            Threads.@spawn let r = $ch[1], i = $ch[2]
                # @info "Chunk $(i), thread $(threadid()), $(Threads.threadpool(threadid())): Spawned $(time() - start)"                
                push!(s, work(r))
                # @info "$(i): Finished $(time() - start)"
            end
        end
    end
    # @info "Finished $(time() - start)"
    # @show s
end
end
using Main.mwe_tasks_2;
mwe_tasks_2.test()
ts = []
for n = 1:10
    push!(ts, @elapsed mwe_tasks_2.test())
end
@show extrema(ts)
