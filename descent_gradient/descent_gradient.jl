#
#   Descent Gradient Method
#
#  x       = Initial guess
#  f       = objective function
#  g       = Gradient of objective function
#  epsilon = tolerance of gradient norm
#  maxiter = maximum of iterations



function descentgradient(x,f,g,epsilon,maxiter,linesearch)
    error = 0
    iter = 0
    stp = -1.0
    while true
        gfx = g(x)
        ngfx = norm(gfx)

        fx = f(x)

        println("$iter $ngfx  $fx  $stp")

        # Convergence checking
        if ngfx < epsilon
            return(x,ngfx,iter,error)
        end

        iter = iter + 1
        # Maximum tteration numbers checking
        if iter > maxiter
            error = 1
            return(x,ngfx,iter,error)
        end

        # Linesearch
        stpmin = 1.e-10
        (stp,x,error) =  linesearch(x,f,gfx,stpmin)

    end


end
