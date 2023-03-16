using Revise; using Pkg; Pkg.activate("."); Pkg.instantiate();
using LinearAlgebra
@show Threads.nthreads()                                                        
@show LinearAlgebra.BLAS.get_num_threads()
include(joinpath(pwd(), "steady_state/3-d/Poisson_examples.jl"))
using .Main.Poisson_examples;
Main.Poisson_examples.Poisson_FE_H20_parass_tasks_example()
Main.Poisson_examples.Poisson_FE_H20_parass_threads_example()

