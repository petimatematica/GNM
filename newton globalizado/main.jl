
using NLPModels, LinearAlgebra, CUTEst

include("gradient.jl")
include("linesearch.jl")
include("newton.jl")
include("hybrid.jl")


    maxiter = 100
    delta = 1.e-2 #troca da direção
    epsilon = 1.e-6
    gamma = 1.e-4
    linesearchG = armijo
    linesearchN = wolfe
    
    
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
