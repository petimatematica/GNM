using NLPModels, LinearAlgebra, CUTEst
include("hybrid.jl")
include("linesearch.jl")

maxiter = 100
delta = 1.e-2 # Troca da direção
epsilon = 1.e-6
gamma = 1.e-4
linesearchG = armijo
linesearchN = wolfe

nlp = CUTEstModel("ROSENBR")
x0 = [1.0, -1.0]

eco_filename = "resultados.txt"

open(eco_filename, "w") do file
    x, error = hybrid()
    println(file, "Para o ponto inicial $x0, o resultado de hybrid foi x = $x e error = $error\n")
end

println("Resultados gravados em $eco_filename")

