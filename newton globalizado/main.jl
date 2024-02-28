#
# Newton method
#

using NLPModels, LinearAlgebra, CUTEst
include("gradient.jl")
include("linesearch.jl")
include("newton.jl")

nlp = CUTEstModel("ROSENBR")


# Objectibve function
function f(x)
    return obj(nlp,x)
    #return 100*(x[2]-x[1]^2)^2+(x[1]-1)^2
end

# Gradient of objective function
function gradf(x)
    return grad(nlp,x)
    #G = zeros(2,1)
    #G[1] = -400*x[1]*(x[2]-x[1]^2)+2*(x[1]-1)
    #G[2] = 200*(x[2]-x[1]^2)
    return G
end

# Hessian of objective functionm
function hessf(x)
    return hess(nlp,x)
end


x0 = nlp.meta.x0
#x0 = [-1.2,1]
#x0 = [1,0.99]

maxiter = 1000000
tol = 1.e-6
gamma = 1.e-4

#Solver calling
linesearch = goldstein

sol,error = newton()




println("$sol")

finalize(nlp)
