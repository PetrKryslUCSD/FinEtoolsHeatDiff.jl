module mwe_tasks
using Base.Threads
function work(r)
    s = 0.0
    for j in r
        s = s + exp(-(j - minimum(r))^2)
    end
    s
end
function test()
    nchunks = 5
    N = 100000000
    chunks = [(((i-1)*N+1:i*N), i) for i = 1:nchunks]
    s = Float64[]
    start = time()
    Threads.@sync begin
        for ch in chunks
            @info "$(ch[2]): Started $(time() - start)"
            Threads.@spawn let r = $ch[1], i = $ch[2]
                @info "$(i): Spawned $(time() - start)"
                push!(s, work(r))
                @info "$(i): Finished $(time() - start)"
            end
        end
    end
    @info "Finished $(time() - start)"
    # @show s
end
end
using Main.mwe_tasks;
ts = []
for n = 1:50
    push!(ts, @elapsed mwe_tasks.test())
end
@show extrema(ts)
