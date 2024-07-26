using NLPModels, LinearAlgebra, CUTEst, Printf
include("hybrid_new.jl")
include("linesearch.jl")

output = open("info.dat","w")

#problems = ["ROSENBR"]
#problems = CUTEst.select(min_var=10, max_var=10)
problems = CUTEst.select(contype="unc",min_var=2,max_var=1000)
length(problems)

DELTA_tests = [1.e-4]

epsilon = 1.e-6
maxiter = 50
gamma = 1.e-3
stpmin = 1.e-12

LG = armijo
LN = armijo
LineSearchs = [armijo]

for p in problems
    try
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
                    x,error,iter,N_count,norm_gradf_x,et = hybrid(x_0,f,g,h,epsilon,maxiter,gamma,delta,LG,LN,stpmin)
                    @printf(       "%15s %10.1e  %8d %5d  %20.15e %10s %10s %15.10e %3d\n",p,delta,iter,N_count,norm_gradf_x,string(LG),string(LN),et,error)
                    @printf(output,"%15s %10.1e  %8d %5d  %20.15e %10s %10s %15.10e %3d\n",p,delta,iter,N_count,norm_gradf_x,string(LG),string(LN),et,error)
                end
            end
        end
        finalize(nlp)
        catch 
        finalize(nlp)
    end 
end
