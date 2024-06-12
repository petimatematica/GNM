#
# Hybrid method
#
function hybrid()
    x = copy(x_init)
    error = 0
    iter = 0
    stp = -1.0

    open(eco_filename, "a") do file
        println(file, "x_0 = $x_init\n")
    end

    while true
        f_x = obj(nlp, x) 
        gradf_x = grad(nlp, x)  
        norm_gradf_x = norm(gradf_x, 2)

        if norm_gradf_x > delta
                println(eco_filename, "G: iter = $iter || grad f(x) || = $norm_gradf_x  f(x) = $f_x  alpha = $stp")

            iter += 1

            if iter > maxiter
                error = 1
                write_to_eco(eco_filename, "Número máximo de iterações atingido.")
                return x, error
            end

            if norm_gradf_x < epsilon
                error = 0
                write_to_eco(eco_filename, "Convergência alcançada.")
                return x, error
            end

            d = -gradf_x
            stp, x, error = linesearchG(x,gradf_x,d,f,gamma)

        else
            write_to_eco(eco_filename, "N: iter = $iter || grad f(x) || = $norm_gradf_x  f(x) = $f_x  alpha = $stp")    

            iter += 1

            if iter > maxiter
                error = 1
                write_to_eco(eco_filename, "Número máximo de iterações atingido.")
                return x, error
            end

            if norm_gradf_x < epsilon
                error = 0
                write_to_eco(eco_filename, "Convergência alcançada.")
                return x, error
            end

            H_x = hess(nlp, x)
            d = -H_x \ gradf_x
            stp, x, error = linesearchN(x,gradf_x,d,f,gamma)
        end
    end
end