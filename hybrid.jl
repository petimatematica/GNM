function hybrid(x :: Vector{Float64},fun :: Function ,grad :: Function ,hess :: Function, 
    epsilon :: Float64, maxiter :: Int64, gamma :: Float64, delta :: Float64,
    linesearchG :: Function , linesearchN :: Function ,stpmin :: Float64)

    iter = 0
    dir = " "
    stp = 0.0
    N_count = 0
    t0 = time()

    while true

        fx_k = fun(x)

        gradf_x = grad(x)

        norm_gradf_x = norm(gradf_x,2)

        if isnan(norm_gradf_x)
            return  return x,4,iter,N_count,norm_gradf_x,Inf
        end

        if norm_gradf_x < epsilon
            et = time() - t0
            return x,0,iter,N_count,norm_gradf_x,et
        end

        iter += 1
        if iter > maxiter
            et = time() - t0\
            return x,1,iter,N_count,norm_gradf_x,Inf
        end

        # Direction select
        if norm_gradf_x > delta 
            dir = "G"
            d = -gradf_x
            stp, x, LS_error = linesearchG(x,gradf_x,grad,d,fx_k,gamma,stpmin)
            if LS_error > 0
                return x,2,iter,N_count,norm_gradf_x,Inf
            end

        else
            try
                dir = "N"
                d = - hess(x) \ gradf_x
                N_count += 1
                stp,x,LS_error = linesearchN(x,gradf_x,grad,d,fx_k,gamma,stpmin)
                if LS_error > 0
                    return x,3,iter,N_count,norm_gradf_x,Inf
                end
            catch e
                dir = "G"
                d = -gradf_x
                stp, x, LS_error = linesearchG(x,gradf_x,grad,d,fx_k,gamma,stpmin)
                if LS_error > 0
                    return x,5,iter,N_count,norm_gradf_x,Inf
                end
            end
            
        end



    end


end
