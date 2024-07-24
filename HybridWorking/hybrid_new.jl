function hybrid(x :: Vector{Float64},fun :: Function ,grad :: Function ,hess :: Function, 
                epsilon :: Float64, maxiter :: Int64, gamma :: Float64, delta :: Float64,
                linesearchG :: Function , linesearchN :: Function )

    iter = 0
    dir = " "
    stp = 0.0
    N_count = 0
    t0 = time()

    while true

        fx_k = fun(x)

        gradf_x = grad(x)

        norm_gradf_x = norm(gradf_x,2)

        #@printf("%5d  %20.15e  %20.15e %20.15e %3s \n",iter,norm_gradf_x,fx_k,stp,dir)
        if norm_gradf_x < epsilon
            et = time() - t0
            return x,0,iter,N_count,norm_gradf_x,et
        end

        iter += 1
        if iter > maxiter
            et = time() - t0
            return x,1,iter,N_count,norm_gradf_x,et
        end

        # Direction select
        if norm_gradf_x > delta 
            dir = "G"
            d = -gradf_x
            stp, x, ~ = linesearchG(x,gradf_x,grad,d,fx_k,gamma)
        else
            dir = "N"
            N_count += 1
            d = - hess(x) \ gradf_x
            stp,x,~ = linesearchN(x,gradf_x,grad,d,fx_k,gamma)
        end



    end
    

end

