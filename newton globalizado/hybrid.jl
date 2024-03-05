#
# Hybrid method
#
function hybrid()
    x = copy(x0)
    error = 0

    # Iteration number setting
    iter = 0
    stp = -1.0

    while true
        # Gradient of f at x
        gradf_x = gradf(x)

        # Norm of gradient of f at x
        norm_gradf_x = norm(gradf_x, 2)

        if norm_gradf_x > 0.25
        
            # Print info
            f_x = f(x)
            println("G: iter = $iter || grad f(x) || = $norm_gradf_x  f(x) = $f_x  alpha = $stp")

            # Convergence checking
            if norm_gradf_x < tol
            error = 0
            return (x, error)
            end

            iter = iter + 1

            # Maximum iteration numbers checking
            if iter > maxiter
            error = 1
            return x, error
            end

            d = -gradf_x
            (stp,x,error) =  linesearchG(x,gradf_x,d,f_x,gamma)

            # Switch to newton method
        
            else
            
            # Print info
            f_x = f(x)
            println("N: iter = $iter || grad f(x) || = $norm_gradf_x  f(x) = $f_x  alpha = $stp")    
           
            # Calculates hessian of f at x
            H_x = hessf(x)
            
            # Convergence checking
            if norm_gradf_x < tol
                error = 0
                return (x, error)
            end
    
                iter = iter + 1
    
            # Maximum iteration numbers checking
            if iter > maxiter
                error = 1
                return x, error
            end
            # Solving Newton equation
            d = -H_x \ gradf_x

            # Linesearch
            stp,x,erro = linesearchN(x,gradf_x,d,f_x,gamma)
        
        end
    end
end