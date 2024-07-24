using NLPModels, LinearAlgebra, CUTEst, Printf
include("hybrid_new.jl")
include("linesearch.jl")

output = open("info.dat","w")

#problems = ["ROSENBR"]
problems = CUTEst.select(min_var=2, max_var=2)
length(problems)

DELTA_tests = [1.0,1.e-2,1.e-4]

epsilon = 1.e-8
maxiter = 10000
gamma = 1.e-3

# LG = armijo
# LN = armijo
LineSearchs = [armijo,goldstein,wolfe]

for p in problems

    nlp = CUTEstModel(p)

    global function f(x)
        return obj(nlp,x)
    end

    global function g(x)
        return grad(nlp,x)
    end

    global function h(x)
        return hess(nlp,x)
    end

    x_0 = nlp.meta.x0
    for LG in LineSearchs
        for LN in LineSearchs

            for delta in DELTA_tests
                x,error,iter,N_count,norm_gradf_x,et = hybrid(x_0,f,g,h,epsilon,maxiter,gamma,delta,LG,LN)
                @printf(       "%9s %10.1e  %5d %5d  %20.15e %10s %10s %15.10e %3d\n",p,delta,iter,N_count,norm_gradf_x,string(LG),string(LN),et,error)
                @printf(output,"%9s %10.1e  %5d %5d  %20.15e %10s %10s %15.10e %3d\n",p,delta,iter,N_count,norm_gradf_x,string(LG),string(LN),et,error)
            end
        end
    end
    
    finalize(nlp)
end


close(output)
