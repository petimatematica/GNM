#
# Newton method
#

using CUTEst, NLPModels, LinearAlgebra
include("gradient.jl")


nlp = CUTEstModel("ROSENBR")


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



maxiter = 1000
tol = 1.e-8
gamma = 1.e-4

#Solver calling
linesearch = armijo

sol,error = gradient()

println("$sol")

finalize(nlp)
