#
# Newton method
#

using CUTEst, NLPModels, LinearAlgebra
include("gradient.jl")
include("linesearch.jl")

nlp = CUTEstModel("HUMPS")


# Objectibve function
function f(x)
    return obj(nlp,x)
end

# Gradient of objective function
function gradf(x)
    return grad(nlp,x)
end

# Hessian of objective functionm
function hessf(x)
    return hess(nlp,x)
end


x0 = nlp.meta.x0

#x0 = [1,0.99]

maxiter = 100
tol = 1.e-6
gamma = 1.e-4

#Solver calling
linesearch = goldstein

sol,error = gradient()

println("$sol")

finalize(nlp)
