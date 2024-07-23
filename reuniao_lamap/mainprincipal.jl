using NLPModels, LinearAlgebra, CUTEst
include("circulo.jl")
include("hybrid.jl")
include("linesearch.jl")

C = [1.0, 2.0]
r = 1
np = 10

points = points_on_circule(C, r, np)


maxiter = 100
delta = 1.e-2 #troca da direção
epsilon = 1.e-6
gamma = 1.e-4
linesearchG = armijo
linesearchN = wolfe

nlp = CUTEstModel("ROSENBR")

function test(points, eco_filename)
    open(eco_filename, "w") do file
        for point in eachrow(points)
            global x_init = point
            x, error = hybrid()
            println(file, "Para o ponto inicial $point, o resultado de hybrid foi x = $x e error = $error\n")
        end
    end
end

eco_filename = "resultados.txt"

test(points, eco_filename)

println("Resultados gravados em $eco_filename")

