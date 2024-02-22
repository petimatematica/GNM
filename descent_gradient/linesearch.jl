#
# Linesearch library
#

#
#  Armijo
#
#
function armijo(x,f,gfx,stpmin)
    sigma = 1.e-3
    error = 0
    stp = 1.0
    fx = f(x)
    gtg = sigma * gfx' * gfx
    while true
        q = x - stp * gfx
        fq = f(q)
        stptest = fq - fx + stp * gtg
        if stptest > 0.0
            stp = stp / 2.0
           if stp < stpmin
                return(stp,q,error)
            end
        else
           return(stp,q,error)
        end
    end
end



#
#  Goldstein
#
function goldstein(x_k,f,gradf_x, d_k)     
        minstep = 1.e-6
        eta1 = 0.25
        eta2 = 1 - eta1   
        gtd = dot(gradf_x,d_k)
        alpha = 1.0;
        fx_k = f(x_k)
        while true
            x_kp1 = x + alpha * d_k
            fx_kp1 = f(x_kp1)

            #Inequalities
            stptestA = ~(fx_kp1 > fx_k + alpha * eta2 * gtd) # Armijo
            stptestB = ~(fx_kp1 < fx_k + alpha * eta1 * gtd)
    
            if stptestA && stptestB  
                return(alpha,x_kp1,0)
            else
                if ~stptestB
                    alpha = eta1 * alpha
                    if alpha < minstep
                        return(alpha,x_kp1,1)
                    end
                else
                    alpha = alpha / eta2
                end
            end
        end
end


#
# Wolfe
#
function wolfe(x,f,g,minstep)   
    eta = 0.25
    eta1 = eta
    eta2 = 1 - eta   
    gTg = gfx' * gfx
    alpha = 1.0;
    fx = f(x)
    while true
        q = x - alpha * gfx
        fq = f(q)
        w = dot(gq, -gfx)
        #Inequalities
        stptestA = ~( w < -eta2 * gTg)
        stptestB = ~(fq - fx > -alpha * eta1 * gTg)

        if stptestA && stptestB  
            return(alpha,q,0)
        else
            if ~stptestB
                alpha = eta1 * alpha
                if alpha < minstep
                    return(alpha,q,1)
                end
            else
                alpha = alpha / eta2
            end
        end
    end
end