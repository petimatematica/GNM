#
#  Main for Descent Gradient Method
#

include("testfunctions.jl")
include("descent_gradient.jl")
include("linesearch.jl")


using LinearAlgebra

# Dimension of problem
n = 2

# Setting objective function 
f(x) = rosenbrock2(x)

# Setting gradient of objective function
g(x) = grad_rosenbrock2(x)

x0 = [2.0;-2.0]
epsilon = 1.e-8
maxiter = 1000
(x,ngfx,iter,error) = descentgradient(x0,f,g,epsilon,maxiter)



