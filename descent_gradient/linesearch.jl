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
#   Goldstein teste
#
function goldstein_test(x,f,gfx,stpmin)

    error = 0
    sigma1 = 1.e-3
    sigma2 = 1.0 - sigma1
    theta1 = 1.5
    theta2 = 0.5

    fx = f(x)
    
    gf = gfx' * gfx
    gf1 = sigma1 * gf
    gf2 = sigma2 * gf

    alpha  =  2.0
    alphaA = 0.0
    alphaV = 1.e-16

    while true
        println("alpha = $alpha")        
        q = x - alpha * gfx
        fq = f(q)

        stptest1 = fq - fx - gf2
        stptest2 = fq - fx - gf1

        flag1 = ~(stptest1 > 0.0) # flag1 = True Armijo OK
        flag2 = ~(stptest2 < 0.0) # flag2 = True Outra ok

        if flag1 && flag2
            return alpha,q,error
        else
            if ~flag1
                alphaA = alpha
            else
                alphaV = alpha
            end
        end
        if alphaA < 1.e-18
            alpha = theta1 * alpha
        else
            alpha = (1.0 - theta2)*alphaV + theta2 * alphaA 
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
