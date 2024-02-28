#
# Armijo 
#
function armijo(x_k,gradf_x,d_k,fx_k,gamma)
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
        end

    end



end
#
# Goldstein
#
function goldstein(x_k,gradf_x,d_k,fx_k,gamma)     
    minstep = 1.e-6
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
                if alpha < minstep
                    return(alpha,x_kp1,1)
                end
            else
                alpha = alpha * eta1
            end
        end
    end
end
