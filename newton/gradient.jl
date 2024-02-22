#
# Gradient method
#
function gradient()
    x = copy(x0)
    error = 0

    # Iteration number seting
    iter = 0
    stp = -1.0
    while true

        # Gradient of f at x
        gradf_x = gradf(x)

        # Norm of gradient of f at x
        norm_gradf_x = norm(gradf_x,2)

        fx = f(x)
        println("$iter $ngfx  $fx  $stp")

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
        (stp,x,error) =  linesearch(x,gradf_x,f_x,gamma)

    end

end

