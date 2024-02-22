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

#
# Armijo Linesearch
#
function armijo(x_k,gradf_x,d_k,fx_k,gamma)
    
    alpha = 1.0
    gtd = gamma * dot(gradf_x,d_k)

    while true
        x_kp1 = x_k + alpha * d_k
        fx_kp1 = f(x_kp1)

        alpha_test = fx_kp1 > fx_k + alpha * gtd

        if ~alpha_test
            return alpha
        else
            alpha = alpha / 2.0
        end

    end



end
