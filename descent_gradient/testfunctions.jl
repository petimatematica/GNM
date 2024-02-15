#
#   Rosenbrock function on R^2
#
function rosenbrock2(x)
    return (x[1]^2 - x[2])^2 + (x[1] - 1.0)^2
end
#
#  Gradient of Rosenbrock function
#
function grad_rosenbrock2(x)
    G1 = 4*x[1] *(x[1]^2 - x[2]) + 2 * (x[1] - 1.0)
    G2 = - 2 * (x[1]^2 - x[2])
    return [G1;G2]
end