#
# Newton method
#
function newton()

    x = copy(x0)

    # Iteration number seting
    iter = 0
    alpha = -1.0

    while true
       
        # Gradient of f at x
        gradf_x = gradf(x)

        # Norm of gradient of f at x
        norm_gradf_x = norm(gradf_x,2)

        # Print info
        f_x = f(x)
        println("iter = $iter || grad f(x) || = $norm_gradf_x  f(x) = $f_x  alpha = $alpha")

        # Detect whether the tolerance was achieved
        if norm_gradf_x < tol
            return x,0
        end
        
        # Check whether the maximum of iterations was achieved
        iter += 1
        if iter > maxiter
            return x,1
        end

        # Calculates hessian of f at x
        H_x = hessf(x)

        # Solving Newton equation
        d = -H_x \ gradf_x

        # Linesearch
        alpha = linesearch(x,gradf_x,d,f_x,gamma)

        x = x + alpha * d
        
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
