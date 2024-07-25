#
# Armijo 
#
function armijo(x_k,gradf_x,g,d_k,fx_k,gamma,stpmin)
    alpha = 1.0
    gtd = gamma * dot(gradf_x,d_k)

    while true
        x_kp1 = x_k + alpha * d_k
        fx_kp1 = f(x_kp1)

        alpha_test = fx_kp1 > fx_k + alpha * gtd

        if ~alpha_test
            return alpha,x_kp1,0
        else
            alpha = alpha / 2.0
            if alpha < stpmin
                return alpha,x_kp1,1
            end
        end

    end



end
#
# Goldstein
#
function goldstein(x_k,gradf_x,g,d_k,fx_k,gamma,minstep)     
    #minstep = 1.e-6
    eta1 = 0.25
    eta2 = 1 - eta1   
    gtd = dot(gradf_x,d_k)
    alpha = 1.0;
   # fx_k = f(x_k)
    while true
        x_kp1 = x_k + alpha * d_k
        fx_kp1 = f(x_kp1)

        #Inequalities
        stptestA = ~(fx_kp1 > fx_k + alpha * eta1 * gtd) # Armijo
        stptestB = ~(fx_kp1 < fx_k + alpha * eta2 * gtd)

        if stptestA && stptestB  
            return(alpha,x_kp1,0)
        else
            if ~stptestB
                alpha =  alpha / eta2
            else
                alpha = alpha * eta1
                if alpha < minstep
                    return(alpha,x_kp1,1)
                end
            end
        end
    end
end
#
# Wolfe
#
function wolfe(x_k,gradf_x,g,d_k,fx_k,gamma,minstep)     
    #minstep = 1.e-6
    eta1 = 0.25
    eta2 = 1 - eta1   
    gtd = dot(gradf_x,d_k)
    alpha = 1.0
    count_test = 0
    while true
        count_test += 1
        if count_test > 10000
            return(alpha,x_k,2)
        end
        #println("$alpha")
        x_kp1 = x_k + alpha * d_k
        fx_kp1 = f(x_kp1)

        #Inequalities
        stptestA = ~(fx_kp1 > fx_k + alpha * eta1 * gtd) # Armijo
        stptestB = ~(dot(g(x_kp1),d_k) < eta2 * gtd)

        if stptestA && stptestB  
            return(alpha,x_kp1,0)
        else
            if ~stptestB
                alpha =  alpha / eta2
                
            else
                alpha = alpha * eta1
                if alpha < minstep
                    return(alpha,x_kp1,1)
                end
            end
        end
    end
end