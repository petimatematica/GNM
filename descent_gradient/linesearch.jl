#
# Linesearch library
#

#
#  Armijo
#
# f(x^k+alpha f'(x^k))<=f(x^k)-sigma alpha||f'(x^k)||^2
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
# Goldstein
# 
#
function goldstein(x,f,gfx,stpmin)   
    error = 0
    beta = 1.e-4
    alpha =  1.0 - beta
    stp = 1.0
    fx = f(x)
    while true
        gtg = gfx' * gfx
        p = x - stp * gfx 
        fp = f(p) 
        steptestone = fp - fx - alpha * gtg * stp #armijo
        steptestwo = fp - fx - beta * gtg * stp 
        if steptestone > 0.0 && steptestwo < 0
            stp = stp / 2.0  
         else 
            if  steptestwo < 0
                stp = stp * alpha
                if stp<stpmin
                   error=1
                   return(stp,p,error) 
                end
            end
            if steptestone > 0.0
                stp = stp / beta
                if stp<stpmin
                    error=1
                    return(stp,p,error) 
                 end 
            end
            return(stp,p,error)
        end    
    end
end


function goldstein_test(x,f,g,minstep)     
        eta = 0.25
        eta1 = eta
        eta2 = 1 - eta   
        gTg = g' * g
        alpha = 1.0;
        fx = f(x)
        while true
            q = x - alpha * g
            fq = f(q)

            #Inequalities
            stptestA = ~(fq - fx < -alpha * eta2 * gTg)
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


#
# Wolfe
#
#
function wolfe(x,f,g,minstep)   
    eta = 0.25
    eta1 = eta
    eta2 = 1 - eta   
    gTg = g' * g
    alpha = 1.0;
    fx = f(x)
    while true
        q = x - alpha * g
        fq = f(q)
        w = dot(grad_rosenbrock2(q), -g)
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