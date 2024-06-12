
using NLPModels, LinearAlgebra, CUTEst

include("gradient.jl")
include("linesearch.jl")
include("newton.jl")
include("hybrid.jl")

function points(r, num_chutes)
    points = []
    for i in 0:num_chutes-1
        theta = 2 * π * i / num_chutes
        x = r * cos(theta)
        y = r * sin(theta)
        push!(points, (x, y))
    end
    return points
end
    maxiter = 100
    delta = 1.e-2 #troca da direção
    epsilon = 1.e-6
    gamma = 1.e-4
    linesearchG = armijo
    linesearchN = wolfe
    
    r = 1.0
    num_chutes = 10
    circle_points = points(r, num_chutes)
    
    eco_filename = "eco_MH.txt"
    open(eco_filename, "w") do file
        println(file, "Número de chutes: $num_chutes\n")
    end
    
    nlp = CUTEstModel("ROSENBR")
    
    for (i, (x, y)) in enumerate(circle_points)
        println("Chute $i: x = $x, y = $y")
        
        x_init = [x, y]
        
        sol, error = hybrid()
        
        open(eco_filename, "a") do file
            println(file, "\nChute $i: x = $x, y = $y")
            if error == 0
                println(file, "Solução: $sol")
            else
                println(file, "Número máximo de iterações atingido.")
            end
        end
    end
    
    finalize(nlp)
