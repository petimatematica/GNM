using Plots
include("circulo.jl")

C = [1.0, 2.0]
r = 1.0
np = 10

T = points_on_circle(C, r, np)

scatter(T)
savefig("circle.png")

println("$T")

