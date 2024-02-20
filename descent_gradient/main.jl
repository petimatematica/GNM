#
#  Main for Descent Gradient Method
#

include("testfunctions.jl")
include("newton_method.jl")
include("linesearch.jl")


using LinearAlgebra

# Dimension of problem
n = 2

# Setting objective function 
f(x) = rosenbrock2(x)

# Setting gradient of objective function
g(x) = grad_rosenbrock2(x)

x0 = [1.0;-1.0]
epsilon = 1.e-8
maxiter = 100
linesearch = Armijo
(x,iter,error) = newtonmethod(x0,f,epsilon,maxiter,linesearch)



