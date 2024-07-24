#
# Gradient method
#
function gradient()
    x = copy(x0)
    error = 0
    iter = 0
    stp = -1.0
    while true
        gradf_x = gradf(x)
        norm_gradf_x = norm(gradf_x,2)

 #       fx = f(x)
 #       println("$iter $norm_gradf_x  $fx  $stp")

        # Print info
        f_x = f(x)
        println("iter = $iter || grad f(x) || = $norm_gradf_x  f(x) = $f_x  alpha = $stp")

        # Convergence checking
        if norm_gradf_x < tol
            error=0
            return(x,error)
        end

        iter = iter + 1

        # Maximum iteration numbers checking
        if iter > maxiter
            error=1
            return(x,error)
        end

        # Linesearch
        d = -gradf_x
        (stp,x,error) =  linesearch(x,gradf_x,d,f_x,gamma)

    end

end

